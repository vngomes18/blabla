# Processador de Arquivos KMZ para Mapas Estáticos

Este script Python converte arquivos KMZ contendo polígonos em mapas estáticos (imagens PNG).

## Requisitos

- Python 3.8 ou superior
- Bibliotecas Python listadas em `requirements.txt`

## Instalação

1. Clone este repositório ou baixe os arquivos
2. Instale as dependências:
```bash
pip install -r requirements.txt
```

## Como Usar

1. Coloque seu arquivo KMZ no mesmo diretório do script
2. Edite a variável `arquivo_kmz` no script `processar_kmz.py` para apontar para seu arquivo KMZ
3. Execute o script:
```bash
python processar_kmz.py
```

O script irá gerar um arquivo `mapa_poligonos.png` contendo o mapa com os polígonos.

## Funcionalidades

- Extrai polígonos de arquivos KMZ
- Converte coordenadas para o formato adequado
- Gera um mapa estático com:
  - Polígonos coloridos
  - Legenda com nomes das áreas
  - Elementos do mapa (costas, fronteiras, etc.)
  - Ajuste automático da visualização

## Personalização

Você pode ajustar os seguintes parâmetros no script:
- Tamanho da figura (`figsize`)
- Cores dos polígonos
- Resolução da imagem de saída (`dpi`)
- Elementos do mapa exibidos 