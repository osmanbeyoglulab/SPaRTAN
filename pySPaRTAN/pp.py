def normalize_column(A, T=0):
    """ perform l2 normalization column-wize of given matrix

    Parameters:
        A : the matrix that works on
        T : switch of column-wize and row-wize.
            T=0: column-wize
            T=1: row-wize
    """

    if (T == 0):
        return np.divide(A, np.sqrt(np.sum(A**2, 0)))
    else:
        At = np.transpose(A)
        return np.transpose(np.divide(At, np.sqrt(np.sum(At**2, 0))))

def clr(X):
    if type(X) is pd.DataFrame:
        x=X.to_numpy()
    else:
        x=X.copy()
    x_norm=np.log1p(x/np.exp(np.mean(np.log1p(x), axis=1)).reshape((-1,1)))
    if type(X) is pd.core.frame.DataFrame:
        x_norm=pd.DataFrame(x_norm, columns=X.columns, index=X.index)
    return x_norm