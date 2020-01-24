import torch, pdb
import torch.nn as nn

class BuildEncoder(nn.Module):
    
    '''
    Builds a convolutional neural network (CNN) that compresses input images 
    into a latent vector.
    
    Args:
        input_shape (list/tuple): shape of input images
        latent_dim (int): dimensionality of latent space
        layers (list): contains integer layer sizes
        activation (callable): instantiated activation function
        addl_convs (int): number of additional convolutions per layer
        use_batchnorm (bool): indicates whether to use batchnorm
    
    Inputs:
        batch (tensor): batch of input images
    
    Returns:
        batch (tensor): batch of latent vectors
    '''
    
    def __init__(self, 
                 input_shape,
                 latent_dim,
                 layers, 
                 activation=None, 
                 addl_convs=0,
                 use_batchnorm=True,
                 dropout_rate=0.0):
        
        super().__init__()
        self.input_shape = input_shape # (channels, height, width)
        self.channels = input_shape[0]
        self.latent_dim = latent_dim
        self.layers = layers
        self.activation = activation if activation is not None else nn.ReLU()
        self.addl_convs = addl_convs
        self.use_batchnorm = use_batchnorm
        self.dropout_rate = dropout_rate
        
        # list of operations in the model
        operations = []
        
        # first pre-process image with conv layers
        operations.append(nn.Conv2d(
            in_channels=self.channels,
            out_channels=self.layers[0],
            kernel_size=3,
            stride=1,
            padding=1))
        if self.use_batchnorm:
            operations.append(nn.BatchNorm2d(self.layers[0]))
        operations.append(self.activation)
        operations.append(nn.Conv2d(
            in_channels=self.layers[0],
            out_channels=self.layers[0],
            kernel_size=3,
            stride=1,
            padding=1))
        if self.use_batchnorm:
            operations.append(nn.BatchNorm2d(self.layers[0]))
        operations.append(self.activation)
        self.channels = self.layers[0]
        
        # loop over layers
        for i, layer in enumerate(layers):
            
            # down sample
            operations.append(nn.Conv2d(
                in_channels=self.channels,
                out_channels=layer,
                kernel_size=3,
                stride=2,
                padding=1))
            
            # batch norm
            if self.use_batchnorm:
                operations.append(nn.BatchNorm2d(layer))
                
            # activation
            operations.append(self.activation)
            
            # dropout
            operations.append(nn.Dropout2d(p=dropout_rate))
            
            # loop over additional convolutions
            for j in range(self.addl_convs):
                
                # convolution
                operations.append(nn.Conv2d(
                    in_channels=layer,
                    out_channels=layer,
                    kernel_size=3,
                    stride=1,
                    padding=1))
            
                # batch norm
                if self.use_batchnorm:
                    operations.append(nn.BatchNorm2d(layer))
                    
                # activation
                operations.append(self.activation)
                
                # dropout
                operations.append(nn.Dropout2d(p=dropout_rate))
                
             # update book keeping
            self.channels = layer
            
        # output linear layer for latent space
        final_shape = int(self.input_shape[1] / 2**len(self.layers))
        self.output_linear = nn.Linear(
            in_features=layer * final_shape**2,
            out_features=self.latent_dim)
                    
        # convert list to sequential model
        self.model = nn.Sequential(*operations)
        
    def forward(self, batch):
        
        # run the fully convolutional model
        batch = self.model(batch)
        
        # flatten and linear
        batch = self.output_linear(batch.view(len(batch), -1))
        
        return batch
        