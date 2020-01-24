import torch, pdb
import torch.nn as nn
import numpy as np

class BuildDecoder(nn.Module):
   
    '''
    Builds a convolutional neural network (CNN) that upscales a latent vector 
    into an image.
   
    Args:
        latent_dim (int): dimensionality of latent space
        init_shape (list/tuple): contains initial image shape
        layers (list): contains integer layer sizes
        output_channels (int): number of output channels
        activation (callable): instantiated activation function
        addl_convs (int):  number of additional convolutions per layer
        use_batchnorm (bool): indicates whether to use batchnorm
   
    Inputs:
        batch (tensor): input batch of latent vectors
   
    Returns:
        batch (tensor): output batch of image reconstructions
    '''
   
    def __init__(self,
                 latent_dim,
                 init_shape,
                 layers,
                 output_channels,
                 activation=None,
                 addl_convs=0,
                 use_batchnorm=True,
                 dropout_rate=0.0):
       
        super().__init__()
        self.latent_dim = latent_dim
        self.init_shape = init_shape # (channels, height, width)
        self.channels = init_shape[0]
        self.layers = layers
        self.output_channels = output_channels
        self.activation = activation if activation is not None else nn.ReLU()
        self.addl_convs = addl_convs
        self.use_batchnorm = use_batchnorm
        self.dropout_rate = dropout_rate
       
        # first convert the latent vector to initial image
        self.input_linear = nn.Linear(
            in_features=self.latent_dim,
            out_features=np.prod(self.init_shape),
            bias=True)
       
        # list of operations in the model
        operations = []
       
        # pre-process initial image with conv layers
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
           
            # upsample
            operations.append(nn.Upsample(
                scale_factor=2, 
                mode='nearest'))
            operations.append(nn.Conv2d(
                in_channels=self.channels,
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
           
        # output convolution
        operations.append(nn.Conv2d(
            in_channels=self.channels,
            out_channels=self.output_channels,
            kernel_size=3,
            stride=1,
            padding=1))
                   
        # convert list to sequential model
        self.model = nn.Sequential(*operations)
       
    def forward(self, batch):
       
        # convert latent vector to initial image
        shape = [len(batch)] + list(self.init_shape)
        batch = self.input_linear(batch).view(shape)
        
        # run the model and squeeze
        batch = self.model(batch)
       
        return batch