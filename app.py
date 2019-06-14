from flask import Flask, render_template, request
import mysql.connector

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('home.html')


@app.route('/database.html', methods=['get'])
def database():
    results = 'hi'
#     connection = mysql.connector.connect(
#         host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
#         user='ophia@hannl-hlo-bioinformatica-mysqlsrv',
#         password='594990',
#         database='ophia')

#     definitions = ['Accessiecode', 'Bit_score', 'Query_coverage',
#                    'Description', 'E_value', 'Header', 'Perc_identity',
#                    'Eiwit_naam', 'Sequentie', 'Taxonomische_Rang']
#     filters = []

#     for defi in definitions:
#         if defi in request.args:
#             filters.append(defi)

#     if not filters:
#         filters = definitions
#     options = ','.join(filters)

#     query = "select " + options + " from eiwitten join resultaten r on " \
#                                   "eiwitten.Eiwit_id = r.Eiwit_id " \
#                                   "join input i " \
#                                   "on r.Sequentie_id = i.Sequentie_id join " \
#                                   "organismen o " \
#                                   "on r.Organisme_id = o.Organisme_id"
#     print(query)

#     cursor = connection.cursor()
#     cursor.execute(query)
#     results = cursor.fetchall()
#     cursor.close()
#     connection.close()
    return render_template('database.html', text=results)


@app.route('/BLAST.html')
def blast():
    return render_template('BLAST.html')


@app.route('/about.html')
def about():
    return render_template('about.html')


if __name__ == '__main__':
    app.run(debug=True)
