import numpy as np

# define the loss function
def loss_function(X, y, theta):
    m = len(y)
    h = np.dot(X, theta)
    loss = np.sum(np.square(h - y)) / (2 * m)
    return loss

# define the gradient function
def gradient(X, y, theta):
    m = len(y)
    h = np.dot(X, theta)
    grad = np.dot(X.T, (h - y)) / m
    return grad

# define the batch gradient descent function
def batch_gradient_descent(X, y, theta, alpha, num_iters):
    m = len(y)
    J_history = np.zeros(num_iters)
    
    for i in range(num_iters):
        grad = gradient(X, y, theta)
        theta = theta - alpha * grad
        J_history[i] = loss_function(X, y, theta)
    
    return theta, J_history

