# fest_params_run.py
from parameters import default_params, drivers_from_params
from model_core import simulate, compute_outputs_along_solution

def main():
    p = default_params()
    drivers = drivers_from_params(p)  # your AnnualCycle/constant-style drivers

    sol = simulate(p, drivers)  # the adapter inside simulate handles it
    series = compute_outputs_along_solution(sol, p, drivers)

    print("OK. Final S/H:", series["S_over_H"][-1])

if __name__ == "__main__":
    main()
    print(default_params())
