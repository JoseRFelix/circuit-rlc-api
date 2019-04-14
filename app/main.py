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

def solve_ODE_equation(L, R, C, t0=None, t1=None, q0=None, i0=None, V=None):
    q = Function("q")
    try:
        E = sympify(V) if V else 0
    except SympifyError:
        return "Wrong input"

    eq = L*q(t).diff(t, 2) + R*q(t).diff(t) + (1/C)*q(t) - E

    condition1 = t0 != None and t1 != None and q0 != None and i0 != None
    condition2 = t0 == None and t1 == None and q0 == None and i0 == None

    if condition1:
        C1, C2, phi = symbols('C1 C2 phi')
        solution_with_constants = dsolve(eq, q(t))

        constants = solve(
            [solution_with_constants.rhs.subs(t, t0)-q0, solution_with_constants.rhs.diff(t).subs(t, t1)-i0])

        solution = solution_with_constants.subs(constants)

        A = sqrt(constants[C1]**2 + constants[C2]**2).evalf(3)
        phi_cos = atan(constants[C2]/constants[C1]).evalf(3)
        phi_sin = atan(constants[C1]/constants[C2]).evalf(3)

        if phi_cos * 57.30 < -90 and sin(phi_cos) > 0 and cos(phi_cos) < 0:
            phi_cos = pi + phi_cos

        if phi_sin * 57.30 < -90 and sin(phi_sin) > 0 and cos(phi_sin) < 0:
            phi_sin = pi + phi_sin

        alternate_solution_cos = trigsimp(dsolve(eq, q(t)).subs(
            {C1: -A*sin(phi), C2: A*cos(phi)})).subs(phi, -phi_cos)

        alternate_solution_sin = trigsimp(dsolve(eq, q(t)).subs(
            {C1: A*cos(phi), C2: A*sin(phi)})).subs(phi, phi_sin)

    elif condition2:
        A, C1, C2, phi, theta = symbols('A C1 C2 phi theta')
        solution_with_constants = dsolve(eq, q(t))

        solution = None

        alternate_solution_cos = None

        alternate_solution_sin = None

    else:
        solution = None
        return solution

    check = checkodesol(eq, solution_with_constants)

    if check[0]:
        return [solution_with_constants, solution, alternate_solution_cos, alternate_solution_sin, None]
    else:
        return [solution_with_constants, solution, alternate_solution_cos, alternate_solution_sin, "No solution"]


class calculate_ODE(Resource):
    def post(self):
        args = parser.parse_args()

        solved_equation = solve_ODE_equation(
            args.L, args.R, args.C, args.t0, args.t1, args.q0, args.i0, args.V)

        if (solved_equation[1] == None):
            solved_equation[1] = latex(None)
        else:
            solved_equation[1] = latex(Eq(solved_equation[1].lhs.diff(t), diff(
                solved_equation[1].rhs, t, 1)))

        result = {'charge_with_constants': latex(solved_equation[0]),
                  'current_with_constants': latex(Eq(solved_equation[0].lhs.diff(t), diff(solved_equation[0].rhs, t, 1))),
                  'charge': latex(solved_equation[1] ),
                  'current':   solved_equation[1],
                  'alternate_solution_cos': latex(solved_equation[2]),
                  'alternate_solution_sin': latex(solved_equation[3])
                  }

        if solved_equation == None:
            return {"Message": "Must send correct initial values"}, 400
        elif solved_equation[4] == "No solution":
            return {**result, "Message": "The DE has no solution"}, 503
        elif solved_equation == "Wrong input":
            return {"Message": "Please input correct function or number"}, 400

        return result, 200


api.add_resource(calculate_ODE, '/calculateODE')

if __name__ == "__main__":
    # Only for debugging while developing
    app.run(host='0.0.0.0', debug=True, port=8000)
