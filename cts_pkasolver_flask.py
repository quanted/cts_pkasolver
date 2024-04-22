from flask import Flask, request, send_file, make_response, send_from_directory, g, Response
import logging
import json
import time

from cts_pkasolver import CTSPkasolver


pkasolver = CTSPkasolver()

app = Flask(__name__)


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("pkasolver")
logger.info("API wrapper for pkasolver.")


@app.route('/pkasolver/test')
def status_check():
    return "pkasolver running.", 200

@app.route('/pkasolver/data')
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
        return {"error": "Missing required pkasolver parater 'smiles'"}, 400
    
    try:
        chart_data, species, pka_list, pka_dict = pkasolver.main(smiles)

        results = {
            "status": True,
            "chart_data": chart_data,
            "species": species,
            "pka_list": pka_list,
            "pka_dict": pka_dict
        }
    except Exception as e:
        logging.error("pkasolver_flask exception: {}".format(e))
        return {"error": "pkasolver internal error"}, 500


    # TODO: Format data for CTS here, or add formatting code on CTS backend?
    # ^^^ Consider what would make sense for implementing future calcs and having
    # to modify the data structure for CTS.

    

    return results, 200



if __name__ == "__main__":
    logging.info("Starting up pkasolver flask app")
    app.run(debug=True, host='0.0.0.0', port=8080)
