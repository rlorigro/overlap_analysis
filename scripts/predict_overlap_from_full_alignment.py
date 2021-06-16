from torch.utils.data.dataset import Dataset
from torch.utils.data import DataLoader
from torch import optim
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

    # print()
    # print(x)
    # print(torch.sigmoid(y_predict).data.numpy().T[:10])
    # print(y[:10])

    y_debug = torch.round(y.cpu().data)
    y_predict_debug = torch.squeeze(y_predict.cpu().data)

    print(y_predict_debug)
    print(y_debug)
    print(torch.round(torch.sigmoid(y_predict_debug)))

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

        y_vectors.append(y.cpu().data.numpy())
        y_predict_vectors.append(torch.sigmoid(y_predict.cpu()).data.numpy())

        batch_index += 1

    y_predict_vector = numpy.concatenate(y_predict_vectors)

    return y_predict_vector


def run(data_loader_train, data_loader_test):
    # n_epochs is the number of times the training set should be iterated over in its entirety
    n_epochs = 2

    # Define the hyperparameters
    learning_rate = 1e-4

    model = AlignmentClassifier()
    model.to(device)
    model.train()

    # Initialize the optimizer with above parameters
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=0.0001)

    # Define the loss function
    loss_fn = nn.BCEWithLogitsLoss()  # Binary cross entropy which does sigmoid 0-1 transform

    # Train and get the resulting loss per iteration
    loss = train(model=model, loader=data_loader_train, optimizer=optimizer, loss_fn=loss_fn, epochs=n_epochs)
    model.eval()

    # Test and get the resulting predicted y values
    y_predict = test(model=model, loader=data_loader_test)

    return loss, y_predict


class AlignmentClassifier(nn.Module):
    def __init__(self):
        # Perform initialization of the pytorch superclass
        super(AlignmentClassifier, self).__init__()

        self.input_size = 4
        self.num_layers = 2
        self.hidden_size = 64
        self.bidirectional = True

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
            bidirectional=self.bidirectional,
            batch_first=True)

        # Define simple linear classifier that operates on flattened RNN output
        in_size = self.num_layers * (1 + self.bidirectional) * self.hidden_size
        out_size = 1

        linear_hidden_size = 128

        self.linear1 = torch.nn.Linear(in_size, linear_hidden_size)
        self.linear2 = torch.nn.Linear(linear_hidden_size, out_size)

    def forward(self,x):
        output,h_n = self.rnn(x)

        # Reorder axes so batch is first ... For some reason still is necessary if "batch_first" is specified for GRU?
        h_n = h_n.permute(1,0,2)

        y = self.linear1(torch.flatten(h_n, start_dim=1))
        y = torch.relu(y)
        y = self.linear2(y)

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
    def __init__(self, directory_path):
        self.directory = os.path.abspath(directory_path)
        self.files = list()
        self.subdirectories = ["",""]
        self.label_indexes = [-1,-1]
        self.load_directories()

    def load_directories(self):
        false_dir = os.path.join(self.directory, "false")
        true_dir = os.path.join(self.directory, "true")

        self.subdirectories[0] = false_dir
        self.subdirectories[1] = true_dir

        if os.path.exists(false_dir):
            self.load_file_paths(false_dir)
        else:
            exit("ERROR: data directory does not contained labeled subdirectory: " + false_dir)

        self.label_indexes[0] = 0
        self.label_indexes[1] = len(self.files)

        if os.path.exists(true_dir):
            self.load_file_paths(true_dir)
        else:
            exit("ERROR: data directory does not contained labeled subdirectory: " + true_dir)

    def load_file_paths(self, subdirectory):
        for path in os.listdir(subdirectory):
            if path.endswith(".csv"):
                self.files.append(path)

    def __len__(self):
        return len(self.files)

    def __getitem__(self, i):
        items = list()

        with open(os.path.join(self.subdirectories[i >= self.label_indexes[1]], self.files[i]),'r') as file:
            for l,line in enumerate(file):
                if l == 0:
                    continue

                if len(line) == 0 or line == '\n':
                    continue

                # kmerId,ordinal0,ordinal1,rlePosition0,rlePosition1,
                item = torch.tensor(list(map(float, line.split(',')[1:-1])), device=device)
                items.append(item)

        x = torch.stack(items, dim=1)
        x = x.permute(1,0)

        # print()
        # print(x[0][1], x[0][0])
        # print(x[:4][:])
        if x[0][1] > x[0][0]:
            p = [1,0,3,2]
            x = x[:,p]
        # print(x[:4][:])

        x -= torch.clone(x[0][:])

        y = torch.zeros([1])
        if i >= self.label_indexes[1]:
            y[0] = 1

        # print(x[:4][:])
        # x = torch.diff(x,dim=0)
        # print(x[:4][:])

        x[:][1:2] /= 100
        x[:][2:] /= 10_000

        return x,y


def main():
    # train_directory = "/home/ryan/data/human/chm13/test"
    # test_directory = "/home/ryan/data/human/chm13/test"

    train_directory = "/home/ryan/data/human/chm13/assembly/subset1/LabeledAlignments"
    test_directory = "/home/ryan/data/human/chm13/assembly/subset2/LabeledAlignments"

    train_dataset = AlignmentDataset(train_directory)
    test_dataset = AlignmentDataset(test_directory)

    n_false = train_dataset.label_indexes[True]
    n_true = len(train_dataset) - n_false

    print(n_false)
    print(n_true)

    false_weight = 1/n_false
    true_weight = 1/n_true

    weights_per_sample = [false_weight]*n_false + [true_weight]*n_true

    print(len(weights_per_sample))

    sampler = torch.utils.data.WeightedRandomSampler(
        weights=weights_per_sample,
        replacement=False,
        num_samples=min(n_false,n_true)*2
    )

    train_dataloader = DataLoader(train_dataset, batch_size=8, collate_fn=pad_collate, sampler=sampler)
    test_dataloader = DataLoader(test_dataset, batch_size=8, collate_fn=pad_collate, shuffle=True)

    run(data_loader_train=train_dataloader, data_loader_test=test_dataloader)


if __name__ == "__main__":
    main()
