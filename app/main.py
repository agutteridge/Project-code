
from flask import Flask, render_template, request, jsonify

from app import controller

# Flask setup
app = Flask(__name__, template_folder='./app/templates',
                      static_folder='./app/static')

@app.route("/")
def index():
    # get another route through which JSON can be sent woooooooooo (maybe same structure for both D3 and GMaps??)
    return render_template('index.html')

@app.route("/data")
def return_data():
    search_term = request.args.get('search_term', '')
    return controller.search_for(search_term)

if __name__ == "__main__":
    app.run(debug=True)
