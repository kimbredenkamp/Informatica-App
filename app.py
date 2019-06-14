

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('home.html')


@app.route('/database.html', methods=['get', 'post'])
def database():
 
    return render_template('database.html')


@app.route('/BLAST.html', methods=['get', 'post'])
def blast():
    param_sequentie = request.args.get('search_word')
    text = ''
    if param_sequentie is not None:
        try:
            result_handle = NCBIWWW.qblast("blastx", "nr", param_sequentie, matrix_name="BLOSUM62", hitlist_size=75,
                                           expect=10, gapcosts="11 2", word_size=6, filter=True)
            blast_records = NCBIXML.read(result_handle)
            for alignment in blast_records.alignments:
                for hsp in alignment.hsps:
                    text += "****Alignment****<br>"
                    text += "sequence:" + str(alignment.title) + "<br>"
                    text += "length:" + str(alignment.length) + "<br>"
                    text += "e value:" + str(hsp.expect) + "<br><br>"
        except TypeError:
            text = "Dit is geen geldige sequentie of accessiecode"

    return render_template('BLAST.html', text=text)


@app.route('/about.html')
def about():
    return render_template('about.html')


if __name__ == '__main__':
    app.run(debug=True)
    
