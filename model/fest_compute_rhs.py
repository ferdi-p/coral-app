# fest_compute_rhs.py
from model_core import compute_rhs_and_internals, make_default_params

def main():
    # params and example state/drivers
    p = make_default_params()
    S, H, jCP, jSG = 1.0, 0.2, 0.2, 0.2
    T, X, Nu = 28.0, 0.5, 0.5

    # base light
    L = 300.0
    d, ints = compute_rhs_and_internals(S, H, jCP, jSG, L, T, X, Nu, p)

    print("\n=== Base call ===")
    print(f"L={L}, T={T}, X={X}, Nu={Nu}")
    print("Derivatives (dS, dH, djCP, djSG):", [f"{v:.6f}" for v in d])
    print("Internals:")
    for k in ["jHG", "jST", "jeL", "jNPQ", "cROS1", "rhoC", "rhoN"]:
        print(f"  {k:6s} = {ints[k]:.6f}")

    # quick sensitivity: increase light
    L2 = 600.0
    d2, ints2 = compute_rhs_and_internals(S, H, jCP, jSG, L2, T, X, Nu, p)

    print("\n=== Light sensitivity ===")
    print(f"jeL:  L={L:>4} -> {ints['jeL']:.6f} | L={L2:>4} -> {ints2['jeL']:.6f}")
    print(f"cROS: L={L:>4} -> {ints['cROS1']:.6f} | L={L2:>4} -> {ints2['cROS1']:.6f}")

if __name__ == "__main__":
    main()

# # run_compute_rhs.py
# from model_core import compute_rhs_and_internals, make_default_params
#
# #make parameters and example inputs
# p = make_default_params()
# S, H, jCP, jSG = 1.0, 0.2, 0.2, 0.2
# L, T, X, Nu = 300.0, 28.0, 0.5, 0.5
#
# # call the function
# d, ints = compute_rhs_and_internals(S, H, jCP, jSG, L, T, X, Nu, p)
#
# # print results
# print("Derivatives (dS, dH, djCP, djSG):", d)
# print("Key internals:")
# for k in ["jHG", "jST", "jeL", "jNPQ", "cROS1", "rhoC", "rhoN"]:
#     print(f"  {k:6s} = {ints[k]:.6f}")
