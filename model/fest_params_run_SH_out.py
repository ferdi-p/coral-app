# fest_params_run.py
from model.parameters import default_params, drivers_from_params
from model.model_core import simulate, compute_outputs_along_solution

def runSH(parameters = default_params()):
    p = parameters
    drivers = drivers_from_params(p)  # your AnnualCycle/constant-style drivers

    sol = simulate(p, drivers)  # the adapter inside simulate handles it
    series = compute_outputs_along_solution(sol, p, drivers)
    #out = series["S_over_H"]
    return series["t"], series["S_over_H"]

print(runSH())

# if __name__ == "__main__":
#     print(runSH())
