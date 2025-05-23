from flask import Flask, render_template, request, send_file, flash, redirect, url_for
import os
from werkzeug.utils import secure_filename
from processar_kmz import processar_arquivo_kmz

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max-limit
app.secret_key = 'sua_chave_secreta_aqui'  # Necessário para flash messages

# Criar pasta de uploads se não existir
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        flash('Nenhum arquivo selecionado', 'error')
        return redirect(url_for('index'))
    
    file = request.files['file']
    if file.filename == '':
        flash('Nenhum arquivo selecionado', 'error')
        return redirect(url_for('index'))
    
    if file and file.filename.endswith('.kmz'):
        try:
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
            
            # Processar o arquivo
            poligonos_validos, poligonos_invalidos = processar_arquivo_kmz(filepath)
            
            # Limpar o arquivo após o processamento
            os.remove(filepath)
            
            # Redirecionar para a página do mapa
            return redirect(url_for('exibir_mapa'))
        except Exception as e:
            flash(f'Erro ao processar arquivo: {str(e)}', 'error')
            return redirect(url_for('index'))
    else:
        flash('Arquivo inválido. Por favor, envie um arquivo .kmz', 'error')
        return redirect(url_for('index'))

@app.route('/mapa')
def exibir_mapa():
    return render_template('mapa.html')

if __name__ == '__main__':
    app.run(debug=True) 