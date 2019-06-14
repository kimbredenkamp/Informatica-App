from flask import Flask, render_template, request

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('home.html')


@app.route('/database.html', methods=['GET', 'POST'])
def database():
    
    return render_template('database.html')


@app.route('/BLAST.html')
def blast():
    

    return render_template('BLAST.html')


@app.route('/about.html')
def about():
    return render_template('about.html')


if __name__ == '__main__':
    app.run(debug=True)
    
