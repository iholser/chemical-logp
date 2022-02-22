from sanic import Sanic
from sanic.response import json

import chem

from typing import Iterable


def _add_cors_headers(response, methods: Iterable[str]) -> None:
    allow_methods = list(set(methods))
    if "OPTIONS" not in allow_methods:
        allow_methods.append("OPTIONS")
    headers = {
        "Access-Control-Allow-Methods": ",".join(allow_methods),
        "Access-Control-Allow-Origin": "*",
        "Access-Control-Allow-Credentials": "true",
        "Access-Control-Allow-Headers": (
            "origin, content-type, accept, "
            "authorization, x-xsrf-token, x-request-id"
        ),
    }
    response.headers.extend(headers)


def add_cors_headers(request, response):
    if request.method != "OPTIONS":
        methods = [method for method in request.route.methods]
        _add_cors_headers(response, methods)


app = Sanic('Holser-logP')
app.register_middleware(add_cors_headers, "response")


@app.route('/')
async def index(request):
    return json({'hello': 'world'})


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
