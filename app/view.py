from flask import Flask, render_template, request, jsonify

from app import controller

# Flask setup
app = Flask(__name__, template_folder='./templates',
                      static_folder='./static')

@app.route("/")
def index():
    return render_template('index.html')

@app.route("/data")
def return_data():
    search_term = request.args.get('search_term', '')
    return controller.search_for(search_term)

def run():
    app.run(debug=True)
