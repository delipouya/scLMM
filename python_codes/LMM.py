import numpy as np
import numpy.linalg as LA


# NOTE:
#   Starting values are all calculated based on var y
#   Later on we can utilize a smarter intialization

# Original MINQE(U,I) n x n formulation
def minqe_unbiased(Kernel, X, y, verbose=True, n_iter=10):
    n = Kernel[0].shape[0]
    p = len(Kernel)
    sig_estimate = np.zeros(p + 1, dtype=np.float64)
    sig_estimate[:] = np.var(y, ddof=1) / (p + 1)
    iteration_estimate = []
    # Convergence
    exit_code = 1

    for i in range(n_iter):
        prev_sig = sig_estimate.copy()

        V = sig_estimate[p] * np.eye(n, dtype=np.float64)  # for I
        for j in range(p):
            V += sig_estimate[j] * Kernel[j]
        # print("\t\tInverting covariance matrix", flush=True)
        V_inv = LA.inv(V)
        R = V_inv - V_inv.dot(X).dot(LA.inv(X.T.dot(V_inv).dot(X))).dot(X.T).dot(V_inv)
        # print("\t\tPopulating Infomation mat and Score vec")
        S = np.zeros((p + 1, p + 1), dtype=np.float64)
        quad = np.zeros(p + 1, dtype=np.float64)
        for j in range(p):
            # print("\t\t\tRow: %d of %d" % (j+1, p), flush=True)
            quad[j] = y.T.dot(R).dot(Kernel[j]).dot(R).dot(y)
            S[j, p] = np.trace(R.dot(Kernel[j]).dot(R))  # last column
            S[p, j] = S[j, p]
            for k in range(j, p):
                S[j, k] = np.trace(R.dot(Kernel[j]).dot(R).dot(Kernel[k]))
                S[k, j] = S[j, k]

        S[p, p] = np.trace(R.dot(R))
        quad[p] = y.T.dot(R).dot(R).dot(y)
        # print("\t\tScore norm: %.6f, abs_max: %.6f" % (LA.norm(quad), np.abs(quad).max()) )
        # print("\t\tSolving system linear system")

        sig_estimate = LA.inv(S).dot(quad)
        iteration_estimate.append(sig_estimate)

        print("Quad:", quad)
        print("Trace (RKRK):", S)

        if verbose:
            print_str = "\t MINQE(U,I) vanilla round " + str(i) + ": "
            for j in range(p + 1):
                print_str += "  {:.4f}".format(sig_estimate[j])
            print(print_str, flush=True)

        diff = np.max(np.abs(sig_estimate - prev_sig))
        if diff < 1e-4:
            exit_code = 0
            if verbose:
                print("\t Estimates converged, exitting - diff: %.6f" % (diff))
            break

    iteration_estimate = np.array(iteration_estimate)
    return sig_estimate, iteration_estimate[0, :], iteration_estimate, exit_code


# Original MINQE(I) n x n formulation
def minqe_biased(Kernel, X, y, verbose=True, n_iter=10):
    n = Kernel[0].shape[0]
    p = len(Kernel)
    sig_estimate = np.zeros(p + 1, dtype=np.float64)
    sig_estimate[:] = np.var(y, ddof=1) / (p + 1)
    iteration_estimate = []
    # Convergence
    exit_code = 1

    for i in range(n_iter):
        prev_sig = sig_estimate.copy()

        V = sig_estimate[p] * np.eye(n, dtype=np.float64)  # for I
        for j in range(p):
            V += sig_estimate[j] * Kernel[j]

        V_inv = LA.inv(V)
        R = V_inv - V_inv.dot(X).dot(LA.inv(X.T.dot(V_inv).dot(X))).dot(X.T).dot(V_inv)

        S = np.zeros((p + 1, p + 1), dtype=np.float64)
        sigma = np.zeros(p + 1, dtype=np.float64)
        for j in range(p):
            sigma[j] = y.T.dot(R).dot(Kernel[j]).dot(R).dot(y)
            S[j, p] = np.trace(V_inv.dot(Kernel[j]).dot(V_inv))  # last column
            S[p, j] = S[j, p]
            for k in range(j, p):
                S[j, k] = np.trace(V_inv.dot(Kernel[j]).dot(V_inv).dot(Kernel[k]))
                S[k, j] = S[j, k]

        S[p, p] = np.trace(V_inv.dot(V_inv))
        sigma[p] = y.T.dot(R).dot(R).dot(y)

        sig_estimate = LA.inv(S).dot(sigma)
        iteration_estimate.append(sig_estimate)

        if verbose:
            print_str = "\t MINQE(I) vanilla round " + str(i) + ": "
            for j in range(p + 1):
                print_str += "  {:.4f}".format(sig_estimate[j])
            print(print_str, flush=True)

        diff = np.max(np.abs(sig_estimate - prev_sig))
        if diff < 1e-4:
            exit_code = 0
            if verbose:
                print("\t Estimates converged, exitting - diff: %.6f" % (diff))
            break

    iteration_estimate = np.array(iteration_estimate)
    return sig_estimate, iteration_estimate[0, :], iteration_estimate, exit_code


# Original EM REML n x n formulation
def reml_em(Kernel, X, y, sig_estimate=None, verbose=True, n_iter=10):
    # Calculating new y
    Q, R = LA.qr(X)
    M = lambda O: O - Q @ (Q.T @ O)
    resid = M(y)
    y_new = (resid - resid.mean()) / resid.std(ddof=1)

    n = Kernel[0].shape[0]
    p = len(Kernel)
    if sig_estimate is None:
        sig_estimate = np.zeros(p + 1, dtype=np.float64)
        sig_estimate[:] = np.var(y_new, ddof=1) / (p + 1)
    iteration_estimate = []
    # Convergence
    exit_code = 1

    for i in range(n_iter):
        prev_sig = sig_estimate.copy()

        V = sig_estimate[p] * np.eye(n, dtype=np.float64)  # for I
        for j in range(p):
            V += sig_estimate[j] * Kernel[j]

        V_inv = LA.inv(V)
        R = V_inv - V_inv.dot(X).dot(LA.inv(X.T.dot(V_inv).dot(X))).dot(X.T).dot(V_inv)

        trace = np.zeros(p + 1, dtype=np.float64)
        quad = np.zeros(p + 1, dtype=np.float64)
        for j in range(p):
            quad[j] = y_new.T.dot(R).dot(Kernel[j]).dot(R).dot(y_new)
            trace[j] = np.trace(R.dot(Kernel[j]))

        quad[p] = y_new.T.dot(R).dot(R).dot(y_new)
        trace[p] = np.trace(R)

        sig_estimate = prev_sig - ((prev_sig ** 2) * (trace - quad)) / n
        iteration_estimate.append(sig_estimate)

        print("Quad:", quad)
        print("Trace (RK):", trace)

        if verbose:
            print_str = "\t EM REML vanilla round " + str(i) + ": "
            for j in range(p + 1):
                print_str += "  {:.4f}".format(sig_estimate[j])
            print(print_str, flush=True)

        diff = np.max(np.abs(sig_estimate - prev_sig))
        if diff < 1e-4:
            exit_code = 0
            if verbose:
                print("\t Estimates converged, exitting - diff: %.6f" % (diff))
            break

    iteration_estimate = np.array(iteration_estimate)
    return sig_estimate, iteration_estimate[0, :], iteration_estimate, exit_code


# Original EM ML n x n formulation
def ml_em(Kernel, X, y, verbose=True, n_iter=10):
    n = Kernel[0].shape[0]
    p = len(Kernel)
    sig_estimate = np.zeros(p + 1, dtype=np.float64)
    sig_estimate[:] = np.var(y, ddof=1) / (p + 1)
    iteration_estimate = []
    # Convergence
    exit_code = 1

    for i in range(n_iter):
        prev_sig = sig_estimate.copy()

        V = sig_estimate[p] * np.eye(n, dtype=np.float64)  # for I
        for j in range(p):
            V += sig_estimate[j] * Kernel[j]

        V_inv = LA.inv(V)
        R = V_inv - V_inv.dot(X).dot(LA.inv(X.T.dot(V_inv).dot(X))).dot(X.T).dot(V_inv)

        trace = np.zeros(p + 1, dtype=np.float64)
        quad = np.zeros(p + 1, dtype=np.float64)
        for j in range(p):
            quad[j] = y.T.dot(R).dot(Kernel[j]).dot(R).dot(y)
            trace[j] = np.trace(V_inv.dot(Kernel[j]))

        quad[p] = y.T.dot(R).dot(R).dot(y)
        trace[p] = np.trace(V_inv)

        sig_estimate = prev_sig - ((prev_sig ** 2) * (trace - quad)) / n
        iteration_estimate.append(sig_estimate)

        if verbose:
            print_str = "\t EM ML vanilla round " + str(i) + ": "
            for j in range(p + 1):
                print_str += "  {:.4f}".format(sig_estimate[j])
            print(print_str, flush=True)

        diff = np.max(np.abs(sig_estimate - prev_sig))
        if diff < 1e-4:
            exit_code = 0
            if verbose:
                print("\t Estimates converged, exitting - diff: %.6f" % (diff))
            break

    iteration_estimate = np.array(iteration_estimate)
    return sig_estimate, iteration_estimate[0, :], iteration_estimate, exit_code


# Original AI REML n x n formulation
def reml_ai(Kernel, X, y, verbose=True, n_iter=10, starting_point=None):
    n = Kernel[0].shape[0]
    p = len(Kernel)
    sig_estimate = np.zeros(p + 1, dtype=np.float64)
    if starting_point is None:
        sig_estimate[:] = np.var(y, ddof=1) / (p + 1)
    else:
        sig_estimate = starting_point.copy()
    iteration_estimate = []
    score_list = []
    information_list = []
    # Convergence
    exit_code = 1

    for i in range(n_iter):
        prev_sig = sig_estimate.copy()

        V = sig_estimate[p] * np.eye(n, dtype=np.float64)  # for I   ## DP: v = sigma_e * I
        print("sig_estimate[-1]:", sig_estimate[p])
        for j in range(p):
            print("sig_estimate[j]:", sig_estimate[j])
            V += sig_estimate[j] * Kernel[j]     ## DP: V = sigma*ZZ' + sigma_e * I
        print("\t\tInverting covariance matrix", flush=True)
        V_inv = LA.inv(V)
        print("\t\tCreating R matrix", flush=True)
        R = V_inv - V_inv.dot(X).dot(LA.inv(X.T.dot(V_inv).dot(X))).dot(X.T).dot(V_inv)

        print("\t\tPopulating Information mat and Score vec")
        Information = np.zeros((p + 1, p + 1), dtype=np.float64)
        score = np.zeros(p + 1, dtype=np.float64)
        for j in range(p):
            print("\t\t\tRow: %d of %d" % (j + 1, p), flush=True)
            score[j] = np.trace(R.dot(Kernel[j])) - y.T.dot(R).dot(Kernel[j]).dot(R).dot(y)
            Information[j, p] = y.T.dot(R).dot(Kernel[j]).dot(R).dot(R).dot(y)
            Information[p, j] = Information[j, p]
            for k in range(j, p):
                Information[j, k] = y.T.dot(R).dot(Kernel[j]).dot(R).dot(Kernel[k]).dot(R).dot(y)
                Information[k, j] = Information[j, k]

        Information[p, p] = y.T.dot(R).dot(R).dot(R).dot(y)
        Information = 0.5 * Information
        score[p] = np.trace(R) - y.T.dot(R).dot(R).dot(y)
        score = -0.5 * score
        print("\t\tScore norm: %.6f, abs_max: %.6f" % (LA.norm(score), np.abs(score).max()))
        score_list.append(score)
        information_list.append(Information)
        print("\t\tSolving system linear system")

        delta = LA.inv(Information).dot(score)
        sig_estimate = sig_estimate + delta
        iteration_estimate.append(sig_estimate)

        print("Information:", Information)
        print("Score:", score)

        if verbose:
            print_str = "\t AI REML vanilla round " + str(i) + ": "
            for j in range(p + 1):
                print_str += "  {:.4f}".format(sig_estimate[j])
            print(print_str, flush=True)

        diff = np.max(np.abs(sig_estimate - prev_sig))
        if diff < 1e-6:
            exit_code = 0
            if verbose:
                print("\t Estimates converged, existing - diff: %.6f" % (diff))
            break

    iteration_estimate = np.array(iteration_estimate)
    score_list = np.array(score_list)
    information_list = np.array(information_list)
    return sig_estimate, iteration_estimate[0, :], iteration_estimate, exit_code, score_list, information_list


# Original AI ML n x n formulation
def ml_ai(Kernel, X, y, verbose=True, n_iter=10):
    n = Kernel[0].shape[0]
    p = len(Kernel)
    sig_estimate = np.zeros(p + 1, dtype=np.float64)
    sig_estimate[:] = np.var(y, ddof=1) / (p + 1)
    iteration_estimate = []
    # Convergence
    exit_code = 1

    for i in range(n_iter):
        prev_sig = sig_estimate.copy()

        V = sig_estimate[p] * np.eye(n, dtype=np.float64)  # for I
        for j in range(p):
            V += sig_estimate[j] * Kernel[j]

        V_inv = LA.inv(V)
        R = V_inv - V_inv.dot(X).dot(LA.inv(X.T.dot(V_inv).dot(X))).dot(X.T).dot(V_inv)

        Information = np.zeros((p + 1, p + 1), dtype=np.float64)
        score = np.zeros(p + 1, dtype=np.float64)
        for j in range(p):
            score[j] = np.trace(V_inv.dot(Kernel[j])) - y.T.dot(R).dot(Kernel[j]).dot(R).dot(y)
            Information[j, p] = y.T.dot(V_inv).dot(Kernel[j]).dot(V_inv).dot(V_inv).dot(y)
            Information[p, j] = Information[j, p]
            for k in range(j, p):
                Information[j, k] = y.T.dot(V_inv).dot(Kernel[j]).dot(V_inv).dot(Kernel[k]).dot(V_inv).dot(y)
                Information[k, j] = Information[j, k]

        Information[p, p] = y.T.dot(V_inv).dot(V_inv).dot(V_inv).dot(y)
        Information = 0.5 * Information
        score[p] = np.trace(V_inv) - y.T.dot(R).dot(R).dot(y)
        score = -0.5 * score

        delta = LA.inv(Information).dot(score)
        sig_estimate = sig_estimate + delta
        iteration_estimate.append(sig_estimate)

        if verbose:
            print_str = "\t AI ML vanilla round " + str(i) + ": "
            for j in range(p + 1):
                print_str += "  {:.4f}".format(sig_estimate[j])
            print(print_str, flush=True)

        diff = np.max(np.abs(sig_estimate - prev_sig))
        if diff < 1e-4:
            exit_code = 0
            if verbose:
                print("\t Estimates converged, exitting - diff: %.6f" % (diff))
            break

    iteration_estimate = np.array(iteration_estimate)
    return sig_estimate, iteration_estimate[0, :], iteration_estimate, exit_code


# Original EI REML n x n formulation
def reml_ei(Kernel, X, y, verbose=True, n_iter=10):
    n = Kernel[0].shape[0]
    p = len(Kernel)
    sig_estimate = np.zeros(p + 1, dtype=np.float64)
    sig_estimate[:] = np.var(y, ddof=1) / (p + 1)
    iteration_estimate = []
    # Convergence
    exit_code = 1

    for i in range(n_iter):
        prev_sig = sig_estimate.copy()

        V = sig_estimate[p] * np.eye(n, dtype=np.float64)  # for I
        for j in range(p):
            V += sig_estimate[j] * Kernel[j]

        V_inv = LA.inv(V)
        R = V_inv - V_inv.dot(X).dot(LA.inv(X.T.dot(V_inv).dot(X))).dot(X.T).dot(V_inv)

        Information = np.zeros((p + 1, p + 1), dtype=np.float64)
        score = np.zeros(p + 1, dtype=np.float64)
        for j in range(p):
            score[j] = np.trace(R.dot(Kernel[j])) - y.T.dot(R).dot(Kernel[j]).dot(R).dot(y)
            Information[j, p] = np.trace(R.dot(Kernel[j]).dot(R))
            Information[p, j] = Information[j, p]
            for k in range(j, p):
                Information[j, k] = np.trace(R.dot(Kernel[j]).dot(R).dot(Kernel[k]))
                Information[k, j] = Information[j, k]

        Information[p, p] = np.trace(R.dot(R))
        Information = 0.5 * Information
        score[p] = np.trace(R) - y.T.dot(R).dot(R).dot(y)
        score = -0.5 * score

        delta = LA.inv(Information).dot(score)
        sig_estimate = sig_estimate + delta
        iteration_estimate.append(sig_estimate)

        print("Information:", Information)
        print("Score:", score)

        if verbose:
            print_str = "\t EI REML vanilla round " + str(i) + ": "
            for j in range(p + 1):
                print_str += "  {:.4f}".format(sig_estimate[j])
            print(print_str, flush=True)

        diff = np.max(np.abs(sig_estimate - prev_sig))
        if diff < 1e-4:
            exit_code = 0
            if verbose:
                print("\t Estimates converged, exitting - diff: %.6f" % (diff))
            break

    iteration_estimate = np.array(iteration_estimate)
    return sig_estimate, iteration_estimate[0, :], iteration_estimate, exit_code


# Original EI ML n x n formulation
def ml_ei(Kernel, X, y, verbose=True, n_iter=10):
    n = Kernel[0].shape[0]
    p = len(Kernel)
    sig_estimate = np.zeros(p + 1, dtype=np.float64)
    sig_estimate[:] = np.var(y, ddof=1) / (p + 1)
    iteration_estimate = []
    # Convergence
    exit_code = 1

    for i in range(n_iter):
        prev_sig = sig_estimate.copy()

        V = sig_estimate[p] * np.eye(n, dtype=np.float64)  # for I
        for j in range(p):
            V += sig_estimate[j] * Kernel[j]

        V_inv = LA.inv(V)
        R = V_inv - V_inv.dot(X).dot(LA.inv(X.T.dot(V_inv).dot(X))).dot(X.T).dot(V_inv)

        Information = np.zeros((p + 1, p + 1), dtype=np.float64)
        score = np.zeros(p + 1, dtype=np.float64)
        for j in range(p):
            score[j] = np.trace(V_inv.dot(Kernel[j])) - y.T.dot(R).dot(Kernel[j]).dot(R).dot(y)
            Information[j, p] = np.trace(V_inv.dot(Kernel[j]).dot(V_inv))
            Information[p, j] = Information[j, p]
            for k in range(j, p):
                Information[j, k] = np.trace(V_inv.dot(Kernel[j]).dot(V_inv).dot(Kernel[k]))
                Information[k, j] = Information[j, k]

        Information[p, p] = np.trace(V_inv.dot(V_inv))
        Information = 0.5 * Information
        score[p] = np.trace(V_inv) - y.T.dot(R).dot(R).dot(y)
        score = -0.5 * score

        delta = LA.inv(Information).dot(score)
        sig_estimate = sig_estimate + delta
        iteration_estimate.append(sig_estimate)

        if verbose:
            print_str = "\t EI ML vanilla round " + str(i) + ": "
            for j in range(p + 1):
                print_str += "  {:.4f}".format(sig_estimate[j])
            print(print_str, flush=True)

        diff = np.max(np.abs(sig_estimate - prev_sig))
        if diff < 1e-4:
            exit_code = 0
            if verbose:
                print("\t Estimates converged, exitting - diff: %.6f" % (diff))
            break

    iteration_estimate = np.array(iteration_estimate)
    return sig_estimate, iteration_estimate[0, :], iteration_estimate, exit_code


print('LMM was imported successfully.')