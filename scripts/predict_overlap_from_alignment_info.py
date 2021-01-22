from Adabound import AdaBoundW, AdaBound
from torch.utils.data.dataset import Dataset
from torch.utils.data import DataLoader
from torch import optim
from matplotlib import pyplot
import torch.nn as nn
import torch
import numpy
import math
import sys

print("USING pytorch VERSION: ", torch.__version__)


def plot_loss(losses, show=True):
    fig = pyplot.gcf()
    fig.set_size_inches(8,6)
    ax = pyplot.axes()
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Loss")
    x_loss = list(range(len(losses)))
    pyplot.plot(x_loss, losses)

    if show:
        pyplot.show()

    pyplot.close()


def load_dataset_as_numpy(csv_path):
    data = list()

    with open(csv_path, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            tokens = line.strip().split(',')

            minAlignedFraction, markerCount, maxDrift, maxSkip, trim, overlap = tokens[2:]

            minAlignedFraction = float(minAlignedFraction)
            markerCount = float(markerCount)
            maxDrift = float(maxDrift)
            maxSkip = float(maxSkip)
            trim = float(trim)
            overlap = int(int(overlap) > 0)

            data.append([minAlignedFraction, markerCount, maxDrift, maxSkip, trim, overlap])

    data = numpy.array(data)

    return data


class ShallowLinear(nn.Module):
    '''
    A simple, general purpose, fully connected network
    '''
    def __init__(self):
        # Perform initialization of the pytorch superclass
        super(ShallowLinear, self).__init__()

        # Define network layer dimensions
        D_in, H1, H2, H3, D_out = [5, 128, 256, 128, 1]    # These numbers correspond to each layer: [input, hidden_1, output]

        # Define layer types
        self.norm0 = torch.nn.BatchNorm1d(D_in)

        self.linear1 = nn.Linear(D_in, H1)
        self.norm1 = torch.nn.BatchNorm1d(H1)

        self.linear2 = nn.Linear(H1, H2)
        self.norm2 = torch.nn.BatchNorm1d(H2)

        self.linear3 = nn.Linear(H2, H3)
        self.norm3 = torch.nn.BatchNorm1d(H3)

        self.linear4 = nn.Linear(H3, D_out)

    def forward(self, x):
        '''
        This method defines the network layering and activation functions
        '''

        # print(x[0])
        x = self.norm0(x)
        # print(x[0])

        x = self.linear1(x) # hidden layer
        x = self.norm1(x)
        x = torch.sigmoid(x)

        x = self.linear2(x) # hidden layer
        x = self.norm2(x)
        x = torch.sigmoid(x)

        x = self.linear3(x) # hidden layer
        x = self.norm3(x)
        x = torch.sigmoid(x)

        x = self.linear4(x) # hidden layer

        return x


def train_batch(model, x, y, optimizer, loss_fn):
    # Run forward calculation
    y_predict = model.forward(x)

    # print()
    # print(x)
    # print(torch.sigmoid(y_predict).data.numpy().T[:10])
    # print(y[:10])

    # Compute loss.
    loss = loss_fn(y_predict.squeeze(), y.float())

    # Before the backward pass, use the optimizer object to zero all of the
    # gradients for the variables it will update (which are the learnable weights
    # of the model)
    optimizer.zero_grad()

    # Backward pass: compute gradient of the loss with respect to model
    # parameters
    loss.backward()

    # Calling the step function on an Optimizer makes an update to its
    # parameters
    optimizer.step()

    return loss.data.item()


def train(model, loader, optimizer, loss_fn, epochs=5):
    losses = list()

    batch_index = 0
    for e in range(epochs):
        for x, y in loader:
            loss = train_batch(model=model, x=x, y=y, optimizer=optimizer, loss_fn=loss_fn)
            print(loss)

            losses.append(loss)

            batch_index += 1

        print("Epoch: ", e+1)
        print("Batches: ", batch_index)

    return losses


def test_batch(model, x, y):
    # run forward calculation
    y_predict = model.forward(x)

    return y, y_predict


def test(model, loader):
    y_vectors = list()
    y_predict_vectors = list()

    batch_index = 0
    for x, y in loader:
        y, y_predict = test_batch(model=model, x=x, y=y)

        y_vectors.append(y.data.numpy())
        y_predict_vectors.append(torch.sigmoid(y_predict).data.numpy())

        batch_index += 1

    y_predict_vector = numpy.concatenate(y_predict_vectors)

    return y_predict_vector


def run(data_loader_train, data_loader_test):
    # Batch size is the number of training examples used to calculate each iteration's gradient
    batch_size_train = 32

    # n_epochs is the number of times the training set should be iterated over in its entirety
    n_epochs = 10

    # Define the hyperparameters
    learning_rate = 1e-2
    shallow_model = ShallowLinear()

    # Initialize the optimizer with above parameters
    optimizer = AdaBound(shallow_model.parameters(), lr=learning_rate)

    # Define the loss function
    loss_fn = nn.BCEWithLogitsLoss()  # Binary cross entropy which does sigmoid 0-1 transform

    # Train and get the resulting loss per iteration
    loss = train(model=shallow_model.eval(), loader=data_loader_train, optimizer=optimizer, loss_fn=loss_fn, epochs=n_epochs)

    # Test and get the resulting predicted y values
    y_predict = test(model=shallow_model, loader=data_loader_test)

    return loss, y_predict


def main(csv_path):
    data = load_dataset_as_numpy(csv_path=csv_path)
    numpy.random.shuffle(data)

    n_false = int(numpy.sum(data[:,-1] == 0))
    n_true = int(numpy.sum(data[:,-1] == 1))

    length = data.shape[0]

    train_length = int(round(length*0.6))
    test_length = length - train_length

    print(n_false)
    print(n_true)
    print(float(n_false)/n_true)

    false_weight = 1/n_false
    true_weight = 1/n_true

    weights = [false_weight, true_weight]

    print(false_weight)
    print(true_weight)

    weights_per_sample = [weights[int(round(i))] for i in data[:,-1]]

    print(weights_per_sample[:100])

    x_train = torch.FloatTensor(data[:train_length,:-1])
    x_test = torch.FloatTensor(data[train_length:,:-1])

    y_train = torch.IntTensor(data[:train_length,-1])
    y_test = torch.IntTensor(data[train_length:,-1])

    weights_train = torch.FloatTensor(weights_per_sample[:train_length])
    weights_test = torch.FloatTensor(weights_per_sample[train_length:])

    dataset_train = torch.utils.data.TensorDataset(x_train, y_train)
    dataset_test = torch.utils.data.TensorDataset(x_test, y_test)

    batch_size = 32
    n_samples = n_false - n_false % batch_size
    # n_samples = n_false + n_true

    train_sampler = torch.utils.data.WeightedRandomSampler(
        weights=weights_train,
        replacement=False,
        num_samples=n_samples)

    test_sampler = torch.utils.data.WeightedRandomSampler(
        weights=weights_test,
        replacement=False,
        num_samples=n_samples)

    data_loader_train = DataLoader(dataset=dataset_train, batch_size=batch_size, sampler=train_sampler)
    data_loader_test = DataLoader(dataset=dataset_test, batch_size=batch_size, sampler=test_sampler)

    losses, y_predict = run(data_loader_train=data_loader_train, data_loader_test=data_loader_test)

    print("Final loss:", sum(losses[-100:])/100)
    plot_loss(losses)

    confusion = numpy.zeros([2,2])

    for i in range(y_predict.shape[0]):
        a = int(round(float(y_predict[i][0])))
        b = int(y_test[i].squeeze())

        # print(y_predict[i][0], y_test[i])
        # print(a, b)

        confusion[a,b] += 1

    fig = pyplot.figure()
    axis = pyplot.axes()

    axis.imshow(confusion)

    pyplot.show()
    pyplot.close()


if __name__ == "__main__":

    if len(sys.argv) != 2:
        exit("ERROR: need to provide 1 argument: path of input GFA")

    main(sys.argv[1])
