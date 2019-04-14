from flask import Flask, request
from flask_restful import reqparse, Resource, Api
from sympy import *
from sympy.abc import t
from flask_cors import CORS
from math import ceil
import sys  # development

app = Flask(__name__)
api = Api(app)
CORS(app)

parser = reqparse.RequestParser()
parser.add_argument(
    'L', type=float, required=True, help='Rate to charge for this resource')
parser.add_argument(
    'R', type=float, required=True, help='Rate to charge for this resource')
parser.add_argument(
    'C', type=float, required=True, help='Rate to charge for this resource')
parser.add_argument(
    't0', type=float, help='Rate to charge for this resource')
parser.add_argument(
    't1', type=float, help='Rate to charge for this resource')
parser.add_argument(
    'q0', type=float, help='Rate to charge for this resource')
parser.add_argument(
    'i0', type=float, help='Rate to charge for this resource')
parser.add_argument(
    'V', type=str, help='Rate to charge for this resource')


def round_expr(expr, num_digits):
    return expr.xreplace({n: round(n, num_digits) for n in expr.atoms(Number)})


def solve_ODE_equation(L, R, C, t0=None, t1=None, q0=None, i0=None, V=None):
    q = Function("q")
    try:
        E = sympify(V) if V else 0
    except SympifyError:
        return None

    eq = L*q(t).diff(t, 2) + R*q(t).diff(t) + (1/C)*q(t) - E

    condition1 = t0 != None and t1 != None and q0 != None and i0 != None
    condition2 = t0 == None and t1 == None and q0 == None and i0 == None

    if condition1:        
        solution_with_constants = round_expr(dsolve(eq, q(t)), 3)
        solution_with_constants_diff = Eq(
            solution_with_constants.lhs.diff(t), diff(solution_with_constants.rhs, t, 1))         

        solution = round_expr(dsolve(eq, q(t), ics={q(t0): q0, q(t).diff(t).subs(t, t1): i0}), 3)
        solution_diff = Eq(solution.lhs.diff(
            t), diff(solution.rhs, t, 1))

    elif condition2:
        A, C1, C2, phi, theta = symbols('A C1 C2 phi theta')
        solution_with_constants = round_expr(dsolve(eq, q(t)), 3)
        solution_with_constants_diff = Eq(
            solution_with_constants.lhs.diff(t), diff(solution_with_constants.rhs, t, 1))

        solution = None
        solution_diff = None

    else:
        solution = None
        return solution

    check = checkodesol(eq, solution_with_constants)

    if check[0]:
        return [solution_with_constants, solution_with_constants_diff, solution, solution_diff, None]
    else:
        return [solution_with_constants, solution_with_constants_diff, solution, solution_diff, "No solution"]


class calculate_ODE(Resource):
    def post(self):
        args = parser.parse_args()

        solved_equation = solve_ODE_equation(
            args.L, args.R, args.C, args.t0, args.t1, args.q0, args.i0, args.V)

        if solved_equation == None:
            return {"message": "Must send correct initial values"}, 400

        result = {'charge_with_constants': latex(solved_equation[0]),
                  'current_with_constants': latex(solved_equation[1]),
                  'charge': latex(solved_equation[2]),
                  'current':   latex(solved_equation[3])
                  }

        if solved_equation[4] == "No solution":
            return {**result, "message": "The DE has no solution"}, 200

        return result, 200


api.add_resource(calculate_ODE, '/calculateODE')

if __name__ == "__main__":
    # Only for debugging while developing
    app.run(host='0.0.0.0', debug=True, port=8000)
