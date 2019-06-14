from flask import Flask, render_template, request
import mysql.connector
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('home.html')


@app.route('/database.html', methods=['GET', 'POST'])
def database():
    connection = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        user='ophia@hannl-hlo-bioinformatica-mysqlsrv',
        password='594990',
        database='ophia')

    params = ['Accessiecode', 'Bit_score', 'Query_coverage',
              'Description', 'E_value', 'Header', 'Perc_identity',
              'Eiwit_naam', 'Taxonomische_Rang', 'Sequentie']

    parameters_to_show_list = []
    for param in params:
        if param in request.args:
            parameters_to_show_list.append(param)

    if not parameters_to_show_list:
        parameters_to_show_list = params

    parameters_to_show = ','.join(parameters_to_show_list)

    search_term = request.args.get('search')
    search_db = request.args.get('search_db')
    where = ''
    if search_term is not None:
        search_term = search_term.replace("'", '"')
        if search_db == 'all':
            where = "where Taxonomische_Rang like '%" + search_term + "%' or Eiwit_naam like '%" + search_term + \
                    "%' or Accessiecode like '%" + search_term + "%'"

        elif search_db == 'organism':
            where = "where Taxonomische_Rang like '%" + search_term + "%' or Accessiecode like '%" + search_term + "%'"

        else:
            where = "where Eiwit_naam like '%" + search_term + "%' or Accessiecode like '%" + search_term + "%'"

    cursor = connection.cursor()
    cursor.execute("""select {} from resultaten 
                      natural join eiwitten 
                      natural join organismen 
                      natural join input
                      {}
                      limit 30""".format(parameters_to_show, where))
    results = cursor.fetchall()

    text = '<table> <tr>'
    for param in parameters_to_show_list:
        text += '<th>' + param + '</th>'
    text += '</tr>'

    for result in results:
        text += '<tr>'
        for parameter in result:
            text += '<td>' + str(parameter) + '</td>'
        text += '</tr>'
    text += '</table>'
    cursor.close()
    connection.close()
    return render_template('database.html', text=text)


@app.route('/BLAST.html', methods=['GET', 'POST'])
def blast():
    param_sequentie = request.args.get('search_word')
    text = ''
    if param_sequentie is not None:
        try:
            result_handle = NCBIWWW.qblast("blastx", "nr", param_sequentie, matrix_name="BLOSUM62", hitlist_size=75,
                                           expect=10, gapcosts="11 1", word_size=6, filter=True)
            blast_records = NCBIXML.read(result_handle)
            for alignment in blast_records.alignments:
                for hsp in alignment.hsps:
                    text += "****Alignment****<br>"
                    text += "sequence:" + str(alignment.title) + "<br>"
                    text += "length:" + str(alignment.length) + "<br>"
                    text += "e value:" + str(hsp.expect) + "<br><br>"
        except TypeError:
            text = "Dit is geen geldige sequentie of accessiecode"
        except ValueError:
            text = "Geef een dna of rna sequentie"
    print(text)
    return render_template('BLAST.html', text=text)


@app.route('/about.html')
def about():
    return render_template('about.html')


if __name__ == '__main__':
    app.run(debug=True)
