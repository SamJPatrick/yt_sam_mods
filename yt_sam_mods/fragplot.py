from functools import wraps
import pydot
import numpy as np
import matplotlib.pyplot as plt
from unyt import unyt_quantity


def clear_graph(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        rv = f(*args, **kwargs)
        args[0].graph = None
        return rv
    return newfunc




class TreePlot:
    """
    Make a simple merger tree plot using pydot and graphviz.

    Parameters
    ----------
    tree : merger tree node :class:`~ytree.data_structures.tree_node.TreeNode`
        The merger tree to be plotted.
    dot_kwargs : optional, dict
        A dictionary of keyword arguments to be passed to pydot.Dot.
        Default: None.
    node_function: optional, function
        A function accepting a single argument of a
        :class:`~ytree.data_structures.tree_node.TreeNode` and returning a
        dictionary of keywords to be given to pydot for creating the node
        object on the plot. This can be used to customize the appearance of
        the nodes. See examples below.
        Default: None.
    edge_function: optional, function
        A function accepting two
        :class:`~ytree.data_structures.tree_node.TreeNode` objects and
        returning a dictionary of keywords to be given to pydot for creating
        the edge object on the plot (the lines connecting halos). This can
        be used to customize the appearance of the edges. See examples below.
        Default: None.

    Attributes
    ----------
    size_field : str
        The field to determine the size of each circle.
        Default: 'mass'.
    size_log : bool
        Whether to scale circle sizes based on log of size field.
        Default: True.
    min_mass : float or unyt_quantity
        The minimum halo mass to be included in the plot. If given
        as a float, units are assumed to be Msun.
        Default: None.
    min_mass_ratio : float
        The minimum ratio between a halo's mass and the mass of the
        main halo to be included in the plot.
        Default: None.

    Examples
    --------

    >>> import ytree
    >>> a = ytree.load("tree_0_0_0.dat")
    >>> p = ytree.TreePlot(a[0])
    >>> p.min_mass = 1e6 # Msun
    >>> p.save()

    >>> # customizing nodes
    >>> import ytree
    >>> def my_node(halo):
    ...     label = f"{halo['uid']}"
    ...     my_kwargs = {"label": label, "fontsize": 8, "shape": "square"}
    ...     return my_kwargs
    >>> a = ytree.load("tree_0_0_0.dat")
    >>> p = ytree.TreePlot(a[0], node_function=my_node)
    >>> p.save()

    >>> # customizing edges
    >>> import ytree
    >>> def my_edge(ancestor, descendent):
    ...     if descendent['mass'] < ancestor['mass']:
    ...         color = 'blue'
    ...     else:
    ...         color = 'black'
    ...     my_kwargs = {"color": color, "penwidth": 5}
    ...     return my_kwargs
    >>> a = ytree.load("tree_0_0_0.dat")
    >>> p = ytree.TreePlot(a[0], edge_function=my_edge)
    >>> p.save()

    """

    _min_dot_size = 0.2
    _max_dot_size = 2
    _min_field_size = None
    _max_field_size = None

    #_size_field = 'fragmentation_instability'
    #_color_field = 'max_metallicity'
    _size_field = None
    _color_field = None
    _size_log = True
    _color_log = False
    _node_cmap = 'inferno'
    _min_mass = None
    _min_mass_ratio = None

    def __init__(self, df, dot_kwargs=None,
                 node_function=None, edge_function=None):
    
        self.df = df
        self.tree = df.tree
        self.node_function = node_function
        self.edge_function = edge_function
        self.dot_kwargs = dict()
        if dot_kwargs is None:
            dot_kwargs = {}
        self.dot_kwargs.update(dot_kwargs)
        self.graph = None


    def save(self, filename=None):
        
        if filename is None:
            filename = f"tree_file.pdf"
        if self.graph is None:
            self._plot()

        suffix = filename[filename.rfind(".")+1:]
        func = getattr(self.graph, f"write_{suffix}", None)
        if func is None:
            raise RuntimeError(f"Cannot save to file format: {suffix}.")

        func(filename)
        return filename

    
    def _plot(self):
        
        self.graph = pydot.Dot(graph_type='graph', **self.dot_kwargs)
        self._plot_children(self.tree)

        
    def _plot_children(self, clump):
        
        graph = self.graph
        my_node = self._plot_node(clump)
        if clump.children is None:
            return
        else :
            children = list(clump.children)

        for child in children:
            if ((self.min_mass is not None) and (child[('clump', 'mass')] < self.min_mass)):
                continue
            if ((self.min_mass_ratio is not None) and \
                (child[('clump', 'mass')] / child.parent[('clump', 'mass')] < self.min_mass_ratio)):
                continue

            child_node = self._plot_node(child)

            if self.edge_function is not None:
                edge_kwargs = self.edge_function(child, clump)
            else:
                edge_kwargs = {"penwidth": 5}

            graph.add_edge(pydot.Edge(my_node, child_node, **edge_kwargs))
            self._plot_children(child)

            
    def _plot_node(self, clump):
        
        graph = self.graph
        node_name = f"{clump.clump_id}"
        my_node = graph.get_node(node_name)

        if len(my_node) == 0:
            if self.node_function is not None:
                node_kwargs = self.node_function(clump)
            else:
                node_kwargs = {'style': 'filled', 'label': ''}
                               #'fillcolor': 'red', 'shape': 'circle'}
                               #'fixedsized': 'true'}
            if self.size_field is not None:
                node_kwargs['width'] = self._size_norm(clump)
            if self.color_field is not None:
                node_kwargs['color'] = self._num_to_color(self._color_norm(clump))
            my_node = pydot.Node(node_name, **node_kwargs)
            graph.add_node(my_node)
        else:
            my_node = my_node[0]

        return my_node

    
    def _size_norm(self, clump):
        
        if self._min_field_size is None:
            tdata = self.df.data[('clump', self.size_field)]
            if self.size_log:
                self._min_field_size = tdata[tdata > 0].min()
            else:
                self._min_field_size = tdata.min()
        nmin = self._min_field_size

        if self._max_field_size is None:
            tdata = self.df.data[('clump', self.size_field)]
            self._max_field_size = tdata.max()
        nmax = self._max_field_size

        fval = clump[('clump', self.size_field)]
        if self.size_log:
            val = np.log(fval / nmin) / np.log(nmax / nmin)
        else:
            val = (fval - nmin) / (nmax - nmin)
        val = np.clip(float(val), 0, 1)

        size = val * (self._max_dot_size - self._min_dot_size) + self._min_dot_size
        return size

    
    def _color_norm(self, clump):
        
        tdata = self.df.data[('clump', self.color_field)]
        if self.color_log:
            self._min_field_color = tdata[tdata > 0].min()
        else:
            self._min_field_color = tdata.min()
        nmin = self._min_field_color

        tdata = self.df.data[('clump', self.color_field)]
        if self.color_log:
            self._max_field_color = tdata[tdata > 0].max()
        else :
            self._max_field_color = tdata.max()
        nmax = self._max_field_color

        fval = clump[('clump', self.color_field)]
        if self.color_log:
            val = np.log(fval / nmin) / np.log(nmax / nmin)
        else:
            val = (fval - nmin) / (nmax - nmin)
        val = np.clip(float(val), 0, 1)
        return val

    
    def _num_to_color(self, num):
        
        cmap = plt.get_cmap(self.node_cmap)
        red, green, blue, alpha = np.rint(np.array(cmap(num)) * 255).astype(int)
        hex_str = f'#{red:2x}{green:2x}{blue:2x}'
        return hex_str

 

    @property
    def min_mass(self):
        """
        The minimum halo mass to be included in the plot.
        """
        return self._min_mass

    @min_mass.setter
    @clear_graph
    def min_mass(self, val):
        if not isinstance(val, unyt_quantity):
            val = unyt_quantity(val, 'Msun')
        self._min_mass = val

    @property
    def min_mass_ratio(self):
        """
        The minimum halo mass to main halo mass.
        """
        return self._min_mass_ratio

    @min_mass_ratio.setter
    @clear_graph
    def min_mass_ratio(self, val):
        self._min_mass_ratio = val

    @property
    def size_field(self):
        """
        The field to determine the size of each circle.
        """
        return self._size_field

    @size_field.setter
    @clear_graph
    def size_field(self, val):
        self._size_field = val

    @property
    def color_field(self):
        """
        The field to determine the color of each circle.
        """
        return self._color_field

    @color_field.setter
    @clear_graph
    def color_field(self, val):
        self._color_field = val


    @property
    def size_log(self):
        """
        Whether to scale circle sizes based on log of size field.
        """
        return self._size_log

    @size_log.setter
    @clear_graph
    def size_log(self, val):
        self._size_log = val

    @property
    def color_log(self):
        """
        Whether to color circle sizes based on log of color field.
        """
        return self._color_log

    @color_log.setter
    @clear_graph
    def color_log(self, val):
        self._color_log = val

    @property
    def node_cmap(self):
        """
        Color map for circles.
        """
        return self._node_cmap

    @node_cmap.setter
    @clear_graph
    def node_cmap(self, val):
        self._node_cmap = val


