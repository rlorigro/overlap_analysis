from torch.utils.data.dataset import Dataset
from torch.utils.data import DataLoader
from torch import optim
from matplotlib import pyplot
import torch.nn as nn
import torch
import numpy
import math
import sys
import os


print("USING pytorch VERSION: ", torch.__version__)
if torch.cuda.is_available():
    device = torch.device('cuda')
    print("CUDA available")
else:
    device = torch.device('cpu')
    print("CUDA NOT available")


def train_batch(model, x, y, optimizer, loss_fn):
    # Run forward calculation
    y_predict = model.forward(x)

    # print()
    # print(x)
    # print(torch.sigmoid(y_predict).data.numpy().T[:10])
    # print(y[:10])

    # Compute loss.
    loss = loss_fn(y_predict.squeeze(), y)

    # Before the backward pass, use the optimizer object to zero all of the
    # gradients for the variables it will update (which are the learnable weights
    # of the model)
    optimizer.zero_grad()

    print("training: ", model.training)

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

        y_vectors.append(y.data.cpu().numpy())
        y_predict_vectors.append(torch.sigmoid(y_predict).data.cpu().numpy())

        batch_index += 1

    y_predict_vector = numpy.concatenate(y_predict_vectors)

    return y_predict_vector


def run(data_loader_train, data_loader_test):
    # Initialize model
    shallow_model = AlignmentClassifier().to(device)
    shallow_model.train()

    # n_epochs is the number of times the training set should be iterated over in its entirety
    n_epochs = 10

    # Define the hyperparameters
    learning_rate = 1e-3

    # Initialize the optimizer with above parameters
    optimizer = torch.optim.Adam(shallow_model.parameters(), lr=learning_rate)

    # Define the loss function
    loss_fn = nn.BCEWithLogitsLoss()  # Binary cross entropy which does sigmoid 0-1 transform

    # Train and get the resulting loss per iteration
    loss = train(model=shallow_model, loader=data_loader_train, optimizer=optimizer, loss_fn=loss_fn, epochs=n_epochs)
    shallow_model.eval()

    # Test and get the resulting predicted y values
    y_predict = test(model=shallow_model, loader=data_loader_test)

    print("Done")

    return loss, y_predict


class AlignmentClassifier(nn.Module):
    def __init__(self):
        # Perform initialization of the pytorch superclass
        super(AlignmentClassifier, self).__init__()

        self.input_size = 5
        self.num_layers = 1
        self.hidden_size = 128

        # Define RNN
        # h_n has shape:
        #   (num_layers * num_directions, batch, hidden_size)
        #
        # Hopefully if "batch_first" is specified, then it will be:
        #   (batch, num_layers * num_directions, hidden_size)
        self.rnn = torch.nn.GRU(
            input_size=self.input_size,
            num_layers=self.num_layers,
            hidden_size=self.hidden_size,
            bidirectional=True,
            batch_first=True)

        # Define simple linear classifier that operates on flattened RNN output
        in_size = self.num_layers * 2 * self.hidden_size
        out_size = 1

        self.linear = torch.nn.Linear(in_size, out_size)

    def forward(self,x):
        output,h_n = self.rnn(x)

        # Reorder axes so batch is first ... For some reason still is necessary if "batch_first" is specified for GRU?
        h_n = h_n.permute(1,0,2)

        y = self.linear(torch.flatten(h_n, start_dim=1))

        return y


'''
Take the x,y tensors returned by the dataset and convert them to a 'packed_padded_sequence'
'''
def pad_collate(items):
    lengths = list()
    sorted_indexes = sorted([i for i in range(len(items))], key=lambda x: len(items[x][0]), reverse=True)

    # Create a batch-first shape with remaining dimensions matching the shape of longest seq in batch
    new_shape = [len(items)] + list(items[sorted_indexes[0]][0].shape)
    x = torch.zeros(new_shape, device=device)
    y = torch.zeros([len(items)], device=device)

    for i in sorted_indexes:
        l = items[i][0].shape[0]

        x[i][0:l][:] = items[i][0]
        y[i] = items[i][1]
        lengths.append(l)

    return torch.nn.utils.rnn.pack_padded_sequence(x, lengths=lengths, batch_first=True), y


class AlignmentDataset:
    def __init__(self):
        self.length = 4
        self.data = [
            [0,1,2,3,4],
            [1,2,3,4,5],
            [2,3,4,5,6],
            [3,4,5,6,7],
            [4,5,6,7,8],
        ]

    def __len__(self):
        return self.length

    def __getitem__(self, i):
        items = list()

        for i in range(0,i+1):
            items.append(torch.tensor(list(map(float, self.data[i])), device=device))

        x = torch.stack(items, dim=1)
        x = x.permute(1,0)

        y = torch.zeros([1])
        if i >= 2:
            y[0] = 1

        return x,y


def main():
    train_dataset = AlignmentDataset()
    test_dataset = AlignmentDataset()

    train_dataloader = DataLoader(train_dataset, batch_size=2, collate_fn=pad_collate, shuffle=True)
    test_dataloader = DataLoader(test_dataset, batch_size=2, collate_fn=pad_collate, shuffle=True)

    run(data_loader_train=train_dataloader, data_loader_test=test_dataloader)


if __name__ == "__main__":
    main()
