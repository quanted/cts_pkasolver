from flask import Flask, request, send_file, make_response, send_from_directory, g, Response
import logging
import json

from cts_pkasolver import CTSPkasolver


pkasolver = CTSPkasolver()

app = Flask(__name__)


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("pkasolver")
logger.info("API wrapper for pkasolver.")


@app.route('/pkasolver/test')
def status_check():
    return "pkasolver running.", 200

@app.route('/pkasolver/data/')
def get_data():
    """
    Returns pka data from pkasolver library.
    """
    t0 = time.time()
    args = request.args
    
    
    # TODO: User input data validation!


    if "smiles" in args:
        smiles = args["smiles"]
    else:
        return "Missing required pkasolver parater 'smiles'", 400
    
    if "data_type" in args:
        data_type = args["data_type"]

    # try:
    chart_data, species = pkasolver.main(smiles, data_type)
    # except Exception as e:
    #     logging.error("pkasolver_flask exception: {}".format(e))
    #     return "pkasolver internal error", 500

    # TODO: Format data for CTS here, or add formatting code on CTS backend?
    # ^^^ Consider what would make sense for implementing future calcs and having
    # to modify the data structure for CTS.

    results = {
        "chart_data": chart_data,
        "species": species
    }

    return results, 200



if __name__ == "__main__":
    logging.info("Starting up pkasolver flask app")
    app.run(debug=True, host='0.0.0.0', port=8080)
