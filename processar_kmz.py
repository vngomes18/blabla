import os
import zipfile
import xml.etree.ElementTree as ET
import folium
from folium.plugins import MarkerCluster
import branca.colormap as cm
import numpy as np
import itertools
from shapely.geometry import Point, Polygon, MultiPoint
import pandas as pd
from geopy.distance import geodesic
from rtree import index

def analisar_kml(kml_content):
    """
    Analisa o conteúdo KML usando ElementTree para entender sua estrutura
    """
    print("\nAnalisando estrutura do KML:")
    root = ET.fromstring(kml_content)
    
    # Define o namespace do KML
    ns = {'kml': 'http://www.opengis.net/kml/2.2'}
    
    # Procura por todos os Placemarks
    placemarks = root.findall('.//kml:Placemark', ns)
    print(f"\nEncontrados {len(placemarks)} Placemarks")
    
    for i, placemark in enumerate(placemarks, 1):
        print(f"\nPlacemark {i}:")
        
        # Procura o nome
        name = placemark.find('kml:name', ns)
        if name is not None:
            print(f"Nome: {name.text}")
        
        # Procura por diferentes tipos de geometria
        polygon = placemark.find('.//kml:Polygon', ns)
        if polygon is not None:
            print("Encontrado Polygon")
            coords = polygon.find('.//kml:coordinates', ns)
            if coords is not None:
                print("Coordenadas encontradas")
                # Imprime as primeiras coordenadas para verificação
                coord_text = coords.text.strip()
                print(f"Primeiras coordenadas: {coord_text[:100]}...")
        else:
            print("Nenhum Polygon encontrado neste Placemark")
            
        # Procura por outros elementos para debug
        for elem in placemark:
            print(f"Elemento encontrado: {elem.tag}")

def extrair_kml(arquivo_kmz):
    """Extrai o conteúdo KML do arquivo KMZ."""
    try:
        with zipfile.ZipFile(arquivo_kmz, 'r') as kmz:
            kml_file = kmz.namelist()[0]
            return kmz.read(kml_file).decode('utf-8')
    except Exception as e:
        print(f"Erro ao extrair KML: {e}")
        return None

def processar_kml(kml_content):
    """Processa o conteúdo KML e extrai os polígonos."""
    if not kml_content:
        return [], []
    
    try:
        root = ET.fromstring(kml_content)
        ns = {'kml': 'http://www.opengis.net/kml/2.2'}
        placemarks = root.findall('.//kml:Placemark', ns)
        
        print(f"\nTotal de Placemarks encontrados: {len(placemarks)}")
        
        poligonos = []
        poligonos_invalidos = []
        for placemark in placemarks:
            nome = placemark.find('kml:name', ns)
            nome = nome.text if nome is not None else "Sem nome"
            coords = placemark.find('.//kml:coordinates', ns)
            if coords is not None and coords.text:
                coord_pairs = coords.text.strip().split()
                pontos = []
                for pair in coord_pairs:
                    try:
                        lon, lat, _ = map(float, pair.split(','))
                        pontos.append([lat, lon])  # Folium espera [lat, lon]
                    except ValueError:
                        continue
                if len(pontos) >= 3:
                    lats = [p[0] for p in pontos]
                    lons = [p[1] for p in pontos]
                    if min(lats) >= -90 and max(lats) <= 90 and min(lons) >= -180 and max(lons) <= 180:
                        poligonos.append({
                            'nome': nome,
                            'coordenadas': pontos
                        })
                        print(f"Polígono processado com sucesso: {nome}")
                    else:
                        poligonos_invalidos.append({
                            'nome': nome,
                            'motivo': 'Coordenadas fora dos limites',
                            'coordenadas': pontos
                        })
                        print(f"Polígono {nome} tem coordenadas fora dos limites válidos")
                else:
                    poligonos_invalidos.append({
                        'nome': nome,
                        'motivo': 'Menos de 3 pontos',
                        'coordenadas': pontos
                    })
                    print(f"Polígono {nome} tem menos de 3 pontos")
            else:
                poligonos_invalidos.append({
                    'nome': nome,
                    'motivo': 'Sem coordenadas',
                    'coordenadas': []
                })
                print(f"Polígono {nome} não tem coordenadas")
        return poligonos, poligonos_invalidos
    except Exception as e:
        print(f"Erro ao processar KML: {e}")
        return [], []

def generate_hexagonal_grid(polygon, radius):
    min_x, min_y, max_x, max_y = polygon.bounds
    dx = radius * np.sqrt(3)
    dy = radius * 2
    points = []
    for x, y in itertools.product(np.arange(min_x, max_x, dx), np.arange(min_y, max_y, dy)):
        if polygon.contains(Point(x, y)):
            points.append((x, y))
        if polygon.contains(Point(x + dx / 2, y + dy / 2)):
            points.append((x + dx / 2, y + dy / 2))
    return points

def limpar_circulos_sobrepostos(circles_data):
    """Remove círculos que estão completamente contidos em outros círculos maiores."""
    # Cria DataFrame com os dados dos círculos
    df_circles = pd.DataFrame(circles_data)
    
    # Ordena o DataFrame com base no raio (do maior para o menor)
    df_circles = df_circles.sort_values(by="raio", ascending=False).reset_index(drop=True)
    
    # Cria uma estrutura de índice espacial (R-tree)
    idx = index.Index()
    
    # Adiciona os círculos ao índice R-tree
    for i, row in df_circles.iterrows():
        lat, lon, r = row["latitude"], row["longitude"], row["raio"]
        # Define o bounding box para o círculo
        bbox = (lon - r / 111320, lat - r / 111320, lon + r / 111320, lat + r / 111320)
        idx.insert(i, bbox)
    
    # Lista para armazenar os círculos não contidos
    filtered_circles = []
    
    # Processa os círculos, verificando se estão contidos
    for idx_row, row in df_circles.iterrows():
        lat1, lon1, r1 = row["latitude"], row["longitude"], row["raio"]
        
        # Encontra os círculos no índice que podem estar contidos
        possible_containments = list(idx.intersection((
            lon1 - r1 / 111320,
            lat1 - r1 / 111320,
            lon1 + r1 / 111320,
            lat1 + r1 / 111320
        )))
        
        contained = False
        for i in possible_containments:
            if i != idx_row:  # Evita comparar o círculo com ele mesmo
                row_f = df_circles.iloc[i]
                lat2, lon2, r2 = row_f["latitude"], row_f["longitude"], row_f["raio"]
                # Verifica se o círculo está contido em outro
                distance = geodesic((lat1, lon1), (lat2, lon2)).meters
                if distance + r1 <= r2:
                    contained = True
                    break
        
        if not contained:
            filtered_circles.append(row.to_dict())
    
    return filtered_circles

def criar_mapa_folium(poligonos):
    if not poligonos:
        print("Nenhum polígono válido para exibir.")
        return
    
    # Centraliza o mapa
    all_lats = [lat for p in poligonos for lat, lon in p['coordenadas']]
    all_lons = [lon for p in poligonos for lat, lon in p['coordenadas']]
    center_lat = np.mean(all_lats)
    center_lon = np.mean(all_lons)
    
    # Configuração dos tiles com atribuições
    tiles = {
        'OpenStreetMap': {
            'url': 'https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',
            'attribution': '© OpenStreetMap contributors'
        },
        'Stamen Terrain': {
            'url': 'https://stamen-tiles-{s}.a.ssl.fastly.net/terrain/{z}/{x}/{y}{r}.png',
            'attribution': 'Map tiles by <a href="http://stamen.com">Stamen Design</a>'
        },
        'CartoDB positron': {
            'url': 'https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',
            'attribution': '© OpenStreetMap contributors & © CartoDB'
        }
    }
    
    m = folium.Map(
        location=[center_lat, center_lon],
        zoom_start=12,
        tiles=tiles['OpenStreetMap']['url'],
        attr=tiles['OpenStreetMap']['attribution']
    )
    
    folium.TileLayer(
        tiles['Stamen Terrain']['url'],
        attr=tiles['Stamen Terrain']['attribution'],
        name='Stamen Terrain'
    ).add_to(m)
    
    folium.TileLayer(
        tiles['CartoDB positron']['url'],
        attr=tiles['CartoDB positron']['attribution'],
        name='CartoDB positron'
    ).add_to(m)
    
    marker_cluster = MarkerCluster().add_to(m)
    colormap = cm.LinearColormap(colors=['red', 'yellow', 'green', 'blue', 'purple'], vmin=0, vmax=len(poligonos))
    
    circles_data = []
    for i, poligono in enumerate(poligonos):
        cor = colormap(i)
        coords = poligono['coordenadas']
        polygon = Polygon([(lon, lat) for lat, lon in coords])
        
        folium.Polygon(
            locations=coords,
            popup=f"<b>{poligono['nome']}</b>",
            color='black',
            weight=2,
            fill=True,
            fill_color=cor,
            fill_opacity=0.4,
            tooltip=poligono['nome']
        ).add_to(m)
        
        # Geração dos círculos
        radius_steps = [1000, 750, 500, 250, 100, 50, 25, 10]
        covered_area = MultiPoint().buffer(0)
        
        for step_idx, radius in enumerate(radius_steps):
            uncovered_area = polygon.difference(covered_area)
            if uncovered_area.is_empty:
                break
                
            if uncovered_area.geom_type == "Polygon":
                uncovered_areas = [uncovered_area]
            else:
                uncovered_areas = list(uncovered_area.geoms)
                
            for gap in uncovered_areas:
                new_circle_points = generate_hexagonal_grid(gap, radius / 111000)
                for point in new_circle_points:
                    circles_data.append({
                        "nome": poligono['nome'],
                        "indice": i,
                        "latitude": point[1],
                        "longitude": point[0],
                        "raio": radius
                    })
                covered_area = covered_area.union(MultiPoint(new_circle_points).buffer(radius / 111000))
        
        # Adiciona o marcador central
        coords_array = np.array(coords)
        centro_lat = np.mean(coords_array[:, 0])
        centro_lon = np.mean(coords_array[:, 1])
        folium.Marker(
            location=[centro_lat, centro_lon],
            popup=f"<b>{poligono['nome']}</b>",
            icon=folium.Icon(color='red', icon='info-sign')
        ).add_to(marker_cluster)
    
    # Limpa os círculos sobrepostos
    filtered_circles = limpar_circulos_sobrepostos(circles_data)
    
    # Adiciona os círculos filtrados ao mapa
    for circle in filtered_circles:
        folium.Circle(
            location=(circle["latitude"], circle["longitude"]),
            radius=circle["raio"],
            color=["green", "yellow", "purple", "orange", "pink", "blue", "black", "black"][radius_steps.index(circle["raio"])],
            fill=True,
            fill_opacity=0.1,
            tooltip=f"Nome: {circle['nome']}, Lat: {circle['latitude']}, Lon: {circle['longitude']}, Índice: {circle['indice']}"
        ).add_to(m)
    
    # Cria a pasta static se não existir
    os.makedirs('static', exist_ok=True)
    
    # Salva os dados dos círculos filtrados em um DataFrame
    df_circles = pd.DataFrame(filtered_circles)
    
    # Configura o pandas para usar ponto e vírgula como separador e vírgula como separador decimal
    df_circles.to_csv('static/circulos_dados.csv', 
                     sep=';', 
                     decimal=',', 
                     index=False,
                     float_format='%.6f')  # Formata os números com 6 casas decimais
    
    folium.LayerControl().add_to(m)
    colormap.add_to(m)
    folium.plugins.Fullscreen().add_to(m)
    folium.plugins.MousePosition().add_to(m)
    
    # Adiciona o código JavaScript para receber mensagens
    m.get_root().html.add_child(folium.Element("""
        <script>
            window.addEventListener('message', function(event) {
                if (event.data.type === 'centerMap') {
                    const map = document.querySelector('#map')._leaflet_map;
                    if (map) {
                        const zoom = Math.log2(360 / (event.data.radius / 111320));
                        map.setView([event.data.lat, event.data.lon], Math.min(Math.max(zoom, 4), 18));
                    }
                }
            });
        </script>
    """))
    
    # Salva o mapa na pasta static
    m.save('static/mapa_poligonos.html')
    print("Mapa interativo gerado com sucesso: static/mapa_poligonos.html")
    print(f"Total de círculos gerados: {len(circles_data)}")
    print(f"Total de círculos após limpeza: {len(filtered_circles)}")

def gerar_relatorio(poligonos_invalidos):
    if not poligonos_invalidos:
        print("Nenhum polígono inválido encontrado.")
        return
    with open('relatorio_poligonos_invalidos.txt', 'w', encoding='utf-8') as f:
        for p in poligonos_invalidos:
            f.write(f"Nome: {p['nome']}\nMotivo: {p['motivo']}\nCoordenadas: {p['coordenadas'][:3]}...\n\n")
    print(f"Relatório de polígonos inválidos salvo em relatorio_poligonos_invalidos.txt ({len(poligonos_invalidos)} registros)")

def processar_arquivo_kmz(arquivo_kmz):
    """Processa um arquivo KMZ e gera o mapa e relatório."""
    kml_content = extrair_kml(arquivo_kmz)
    if not kml_content:
        raise Exception("Não foi possível extrair o conteúdo KML do arquivo")
    
    poligonos, poligonos_invalidos = processar_kml(kml_content)
    print(f"\nTotal de polígonos processados com sucesso: {len(poligonos)}")
    print(f"Total de polígonos inválidos: {len(poligonos_invalidos)}")
    
    criar_mapa_folium(poligonos)
    gerar_relatorio(poligonos_invalidos)
    
    return len(poligonos), len(poligonos_invalidos)

def main():
    arquivo_kmz = 'areas.kmz'
    if not os.path.exists(arquivo_kmz):
        print(f"Erro: O arquivo {arquivo_kmz} não foi encontrado.")
        return
    kml_content = extrair_kml(arquivo_kmz)
    if not kml_content:
        return
    poligonos, poligonos_invalidos = processar_kml(kml_content)
    print(f"\nTotal de polígonos processados com sucesso: {len(poligonos)}")
    print(f"Total de polígonos inválidos: {len(poligonos_invalidos)}")
    criar_mapa_folium(poligonos)
    gerar_relatorio(poligonos_invalidos)

if __name__ == "__main__":
    main() 