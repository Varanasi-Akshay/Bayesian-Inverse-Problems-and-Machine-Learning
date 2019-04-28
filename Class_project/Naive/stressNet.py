import torch
from torch import autograd
import torch.nn as nn
import torch.nn.functional as functional
import torch.optim as optim
import sys
from os.path import dirname, realpath
parent = dirname(dirname(realpath(__file__)))
parent+='/sampleGenerator'
sys.path.append(parent)
from sampleGenerator import sampleGenerator
from sklearn.preprocessing import StandardScaler
import IPython

"""
input-strains-output-stress neural network

[1142] loss:     19.508
    rloss: 0.00528005
"""

# Customize the loss func
# https://spandan-madan.github.io/A-Collection-of-important-tasks-in-pytorch/

class stressNet(nn.Module):

    def __init__(self, neurons):
        
        super(stressNet, self).__init__()

        # check inputs for init
        if not isinstance(neurons, list):
            raise TypeError('neorons must be list of integers')
        for i in neurons:
            if not isinstance(i, int):
                raise TypeError('neorons must be list of integers')

        # stores neurons per layer list
        self.neurons = neurons
        self.num_hiddens = len(neurons)

        # an affine operation: y = Wx + b
        setattr(self, 'fc1', nn.Linear(3, self.neurons[0]))

        for i in range(self.num_hiddens-1):
            setattr(
                self, 
                'fc'+str(i+2), 
                nn.Linear(self.neurons[i], self.neurons[i+1])
            )

        setattr(
            self, 
            'fc'+str(self.num_hiddens+1), 
            nn.Linear(self.neurons[-1], 3)
        )

    def forward(self, x):

        for i in range(self.num_hiddens):
            fc = getattr(self, 'fc'+str(i+1))
            x = functional.softplus(fc(x))
        
        fc = getattr(self, 'fc'+str(self.num_hiddens+1))
        x = fc(x)

        return x



def cal_stiff(S, E):
    C = torch.stack(
        [autograd.grad(S[i], E, create_graph=True)[0]\
        for i in range(3)]
    )
    return C


neurons = [10]*2
stress_net = stressNet(neurons)
stress_net.double()

data = sampleGenerator()
e, s = data

scaler_E = StandardScaler()
scaler_S = StandardScaler()

scaler_E.fit(e)
scaler_S.fit(s)
