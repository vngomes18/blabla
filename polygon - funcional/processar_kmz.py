import os
import zipfile
import xml.etree.ElementTree as ET
import folium
from folium.plugins import MarkerCluster
import branca.colormap as cm
import numpy as np

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

def criar_mapa_folium(poligonos):
    if not poligonos:
        print("Nenhum polígono válido para exibir.")
        return
    # Centraliza o mapa
    all_lats = [lat for p in poligonos for lat, lon in p['coordenadas']]
    all_lons = [lon for p in poligonos for lat, lon in p['coordenadas']]
    center_lat = np.mean(all_lats)
    center_lon = np.mean(all_lons)
    m = folium.Map(location=[center_lat, center_lon], zoom_start=12, tiles='OpenStreetMap')
    folium.TileLayer('Stamen Terrain').add_to(m)
    folium.TileLayer('CartoDB positron').add_to(m)
    marker_cluster = MarkerCluster().add_to(m)
    colormap = cm.LinearColormap(colors=['red', 'yellow', 'green', 'blue', 'purple'], vmin=0, vmax=len(poligonos))
    for i, poligono in enumerate(poligonos):
        cor = colormap(i)
        folium.Polygon(
            locations=poligono['coordenadas'],
            popup=f"<b>{poligono['nome']}</b>",
            color='black',
            weight=2,
            fill=True,
            fill_color=cor,
            fill_opacity=0.4,
            tooltip=poligono['nome']
        ).add_to(m)
        coords = np.array(poligono['coordenadas'])
        centro_lat = np.mean(coords[:, 0])
        centro_lon = np.mean(coords[:, 1])
        folium.Marker(
            location=[centro_lat, centro_lon],
            popup=f"<b>{poligono['nome']}</b>",
            icon=folium.Icon(color='red', icon='info-sign')
        ).add_to(marker_cluster)
    folium.LayerControl().add_to(m)
    colormap.add_to(m)
    folium.plugins.Fullscreen().add_to(m)
    folium.plugins.MousePosition().add_to(m)
    m.save('mapa_poligonos.html')
    print("Mapa interativo gerado com sucesso: mapa_poligonos.html")

def gerar_relatorio(poligonos_invalidos):
    if not poligonos_invalidos:
        print("Nenhum polígono inválido encontrado.")
        return
    with open('relatorio_poligonos_invalidos.txt', 'w', encoding='utf-8') as f:
        for p in poligonos_invalidos:
            f.write(f"Nome: {p['nome']}\nMotivo: {p['motivo']}\nCoordenadas: {p['coordenadas'][:3]}...\n\n")
    print(f"Relatório de polígonos inválidos salvo em relatorio_poligonos_invalidos.txt ({len(poligonos_invalidos)} registros)")

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