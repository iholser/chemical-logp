from sanic import Sanic
from sanic.response import json
from sanic_limiter import Limiter, get_remote_address
from typing import Iterable

import chem
app = Sanic('Holser-logP')

limiter = Limiter(app, global_limits=[
                  '200 per minute', '2000 per day'], key_func=get_remote_address)


# def _add_cors_headers(response, methods: Iterable[str]) -> None:
#     allow_methods = list(set(methods))
#     if "OPTIONS" not in allow_methods:
#         allow_methods.append("OPTIONS")
#     headers = {
#         "Access-Control-Allow-Methods": ",".join(allow_methods),
#         "Access-Control-Allow-Origin": "*",
#         "Access-Control-Allow-Credentials": "true",
#         "Access-Control-Allow-Headers": (
#             "origin, content-type, accept, "
#             "authorization, x-xsrf-token, x-request-id"
#         ),
#     }
#     response.headers.extend(headers)


def add_cors_headers(request, response):
    if request.method != "OPTIONS":
        methods = list(
            set(["OPTIONS", *(method for method in request.route.methods)]))
        headers = {
            "Access-Control-Allow-Methods": ",".join(methods),
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Credentials": "true",
            "Access-Control-Allow-Headers": (
                "origin, content-type, accept, "
                "authorization, x-request-id"
            ),
        }
        response.headers.extend(headers)


app.register_middleware(add_cors_headers, "response")


@app.route('/')
async def index(request):
    return json({'hello': 'world'})


@app.route('/smiles', methods=['POST'])
def mol_from_smiles(request):
    try:
        mol = chem.mol_from_smiles(request.body)
        return json({"mol": mol})
    except Exception as e:
        print(e)
        return json({"error": e})


@app.route('/mol', methods=['POST'])
def chem_info(request):
    try:
        # print(request.body)
        mol = chem.chem_3d_from_mol_block(request.body)
        return json(chem.get_chemical_info(mol))
    except Exception as e:
        print(e)
        return json({"error": e})


if __name__ == '__main__':
    app.run(host='0.0.0.0')
