import torch
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
from potentialNet import potential_net, scaler_E, scaler_S, cal_stress
import IPython

"""
input-strains-output-potential neural network

[1142] loss:     19.508
    rloss: 0.00528005
"""

# Customize the loss func
# https://spandan-madan.github.io/A-Collection-of-important-tasks-in-pytorch/


e, s = sampleGenerator()

E = scaler_E.transform(e)
S = scaler_S.transform(s)

E, S = torch.from_numpy(E), torch.from_numpy(S)

E.requires_grad = True
# S.requires_grad = True

criterion = nn.MSELoss()
# optimizer = optim.ASGD(
#     potential_net.parameters(),
#     lr = 2e-3,
#     weight_decay=1e-5,
#     t0 = 1e15,
# )
optimizer = optim.Adam(potential_net.parameters(), lr=1e-3, amsgrad=True)

running_loss = 0.0

num_of_samples = len(E)

running_loss = 0.0
num_of_samples = len(E)
epoch = 0
rloss = 1.0
init_loss = 500.
rloss_old = rloss
rloss_best = rloss
isbest = True

filename='trained_net.pt'
for paras in potential_net.parameters():
    nn.init.normal_(paras)
# potential_net.load_state_dict(torch.load(filename))


while rloss > 0.01 or running_loss < 1000:  # loop over the dataset multiple times
    running_loss = 0.0
    # fix zero point
    x = torch.zeros_like(E[0])
    x.requires_grad = True
    f = potential_net(x)
    loss = criterion(f, torch.zeros_like(f))
    loss.backward()
    for i in range(num_of_samples):
        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        x = E[i]
        f = potential_net(x)
        g = cal_stress(f, x)

        #TODO: grad(outputs) <= outputs
        loss = criterion(g, S[i])
        loss.backward()
        running_loss += loss.item()
        optimizer.step()

    if epoch == 0 and init_loss==0.:
        init_loss = running_loss
    if epoch >= 1:
        rloss = running_loss / init_loss
    
    
    if rloss <= rloss_best:
        rloss_best = rloss
        torch.save(potential_net.state_dict(),filename)

    if epoch==0:
        print('[%d] loss: %10.3f' %(epoch, init_loss))
        print('    rloss: %10.8f' %(rloss))
    if epoch>=1:
        print('[%d] loss: %10.3f' %(epoch, running_loss))
        print('    rloss: %10.8f' %(rloss))
        print('  descent: '+str(rloss < rloss_old))
        print('bestrloss: %10.8f' %(rloss_best))
        print('   isbest: '+str(isbest))
    rloss_old = rloss
    epoch += 1


print('Finished Training')

IPython.embed()
