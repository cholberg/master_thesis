import numpy as np
import torch
from torch import optim
from utils import *

# Get from https://github.com/kimiandj/gsw
class GSW:
    def __init__(self, ftype="linear", degree=2, radius=1000, use_cuda=True):
        self.ftype = ftype
        self.degree = degree
        self.radius = radius
        if torch.cuda.is_available() and use_cuda:
            self.device = torch.device("cuda")
        else:
            self.device = torch.device("cpu")
        self.theta = None  # This is for max-GSW

    def gsw(self, X, Y, theta=None, p=1, L=1000):
        """
        Calculates GSW between two empirical distributions.
        Note that the number of samples is assumed to be equal
        (This is however not necessary and could be easily extended
        for empirical distributions with different number of samples)
        """
        N, dn = X.shape
        M, dm = Y.shape
        assert dn == dm and M == N
        if theta is None:
            theta = self.random_slice(dn, L)

        Xslices = self.get_slice(X, theta)
        Yslices = self.get_slice(Y, theta)

        Xslices_sorted = torch.sort(Xslices, dim=0)[0]
        Yslices_sorted = torch.sort(Yslices, dim=0)[0]
        return (
            (torch.sum(torch.abs(Xslices_sorted - Yslices_sorted) ** p) / N) / L
        ) ** (1 / p)

    def max_gsw(self, X, Y, p=1, iterations=50, lr=1e-4, num_start=30):
        N, dn = X.shape
        M, dm = Y.shape
        device = self.device
        assert dn == dm and M == N
        scores = np.zeros(num_start)
        for i in range(num_start):
            if self.ftype == "linear":
                theta = torch.randn((1, dn), device=device, requires_grad=True)
                theta.data /= torch.sqrt(torch.sum((theta.data) ** 2))
            elif self.ftype == "poly":
                dpoly = self.homopoly(dn, self.degree)
                theta = torch.randn((1, dpoly), device=device, requires_grad=True)
                theta.data /= torch.sqrt(torch.sum((theta.data) ** 2))
            elif self.ftype == "circular":
                theta = torch.randn((1, dn), device=device, requires_grad=True)
                theta.data /= torch.sqrt(torch.sum((theta.data) ** 2))

            optimizer = optim.Adam([theta], lr=lr)
            score = np.zeros((iterations,))
            for j in range(iterations):
                optimizer.zero_grad()
                loss = -self.gsw(
                    X.to(self.device), Y.to(self.device), theta.to(self.device), p, 1
                )
                score[j] = loss.item()
                loss.backward(retain_graph=True)
                optimizer.step()
                theta.data /= torch.sqrt(torch.sum(theta.data ** 2))

            scores[i] = -np.min(score)
        return np.max(scores)

    def gsl2(self, X, Y, theta=None):
        """
        Calculates GSW between two empirical distributions.
        Note that the number of samples is assumed to be equal
        (This is however not necessary and could be easily extended
        for empirical distributions with different number of samples)
        """
        N, dn = X.shape
        M, dm = Y.shape
        assert dn == dm and M == N
        if theta is None:
            theta = self.random_slice(dn)

        Xslices = self.get_slice(X, theta)
        Yslices = self.get_slice(Y, theta)

        Yslices_sorted = torch.sort(Yslices, dim=0)

        return torch.sqrt(torch.sum((Xslices - Yslices) ** 2))

    def get_slice(self, X, theta):
        """Slices samples from distribution X~P_X
        Inputs:
            X:  Nxd matrix of N data samples
            theta: parameters of g (e.g., a d vector in the linear case)
        """
        if self.ftype == "linear":
            return self.linear(X, theta)
        elif self.ftype == "poly":
            return self.poly(X, theta)
        elif self.ftype == "circular":
            return self.circular(X, theta)
        else:
            raise Exception("Defining function not implemented")

    def random_slice(self, dim, L):
        if self.ftype == "linear":
            theta = torch.randn((L, dim))
            theta = torch.stack([th / torch.sqrt((th ** 2).sum()) for th in theta])
        elif self.ftype == "poly":
            dpoly = self.homopoly(dim, self.degree)
            theta = torch.randn((L, dpoly))
            theta = torch.stack([th / torch.sqrt((th ** 2).sum()) for th in theta])
        elif self.ftype == "circular":
            theta = torch.randn((L, dim))
            theta = torch.stack([th / torch.sqrt((th ** 2).sum()) for th in theta])
        return theta.to(self.device)

    def linear(self, X, theta):
        if len(theta.shape) == 1:
            return torch.matmul(X, theta)
        else:
            return torch.matmul(X, theta.t())

    def poly(self, X, theta):
        """The polynomial defining function for generalized Radon transform
        Inputs
        X:  Nxd matrix of N data samples
        theta: Lxd vector that parameterizes for L projections
        degree: degree of the polynomial
        """
        N, d = X.shape
        assert theta.shape[1] == self.homopoly(d, self.degree)
        powers = list(self.get_powers(d, self.degree))
        HX = torch.ones((N, len(powers))).to(self.device)
        for k, power in enumerate(powers):
            for i, p in enumerate(power):
                HX[:, k] *= X[:, i] ** p
        if len(theta.shape) == 1:
            return torch.matmul(HX, theta)
        else:
            return torch.matmul(HX, theta.t())

    def circular(self, X, theta):
        """The circular defining function for generalized Radon transform
        Inputs
        X:  Nxd matrix of N data samples
        theta: Lxd vector that parameterizes for L projections
        """
        N, d = X.shape
        if len(theta.shape) == 1:
            return torch.sqrt(torch.sum((X - self.radius * theta) ** 2, dim=1))
        else:
            return torch.stack(
                [
                    torch.sqrt(torch.sum((X - self.radius * th) ** 2, dim=1))
                    for th in theta
                ],
                1,
            )

    def get_powers(self, dim, degree):
        """
        This function calculates the powers of a homogeneous polynomial
        e.g.

        list(get_powers(dim=2,degree=3))
        [(0, 3), (1, 2), (2, 1), (3, 0)]

        list(get_powers(dim=3,degree=2))
        [(0, 0, 2), (0, 1, 1), (0, 2, 0), (1, 0, 1), (1, 1, 0), (2, 0, 0)]
        """
        if dim == 1:
            yield (degree,)
        else:
            for value in range(degree + 1):
                for permutation in self.get_powers(dim - 1, degree - value):
                    yield (value,) + permutation

    def homopoly(self, dim, degree):
        """
        calculates the number of elements in a homogeneous polynomial
        """
        return len(list(self.get_powers(dim, degree)))
