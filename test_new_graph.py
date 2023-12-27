import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


# create dictionary to swap from idx to idc and vice versa
def create_idc_dictionary(nx_g):
    edge_dict = {}
    for edge_idx, edge_idc in enumerate(nx_g.edges()):
        edge_dict[edge_idx] = tuple(sorted(edge_idc, key=sum))
    return edge_dict


def get_idx_from_idc(edge_dictionary, idc):
    idc = tuple(sorted(idc, key=sum))
    return list(edge_dictionary.values()).index(idc)


def get_idc_from_idx(edge_dictionary, idx):
    return edge_dictionary[idx]


class GraphCreator:
    def __init__(self, m, n, ion_chain_size_vertical, ion_chain_size_horizontal):
        self.m = m
        self.n = n
        self.ion_chain_size_vertical = ion_chain_size_vertical
        self.ion_chain_size_horizontal = ion_chain_size_horizontal
        self.num_edges = self.n // 2
        self.networkx_graph = self.create_graph()

        self.idc_dict = create_idc_dictionary(self.networkx_graph)
        self.pz_edges_idx = [
            get_idx_from_idc(self.idc_dict, edge)
            for edge in self.networkx_graph.edges()
            if nx.get_edge_attributes(self.networkx_graph, "edge_type")[edge] != "trap"
        ]

    def create_graph(self):
        self.m_extended = self.m + (self.ion_chain_size_vertical - 1) * (self.m - 1)
        self.n_extended = self.n + (self.ion_chain_size_horizontal - 1) * (self.n - 1)

        networkx_graph = nx.grid_2d_graph(self.m_extended, self.n_extended)
        self._set_trap_nodes(networkx_graph)
        self._remove_horizontal_edges(networkx_graph)
        self._remove_vertical_edges(networkx_graph)
        self._remove_horizontal_nodes(networkx_graph)
        self._set_junction_nodes(networkx_graph)
        nx.set_edge_attributes(networkx_graph, "trap", "edge_type")
        self._set_processing_zone(networkx_graph)

        return networkx_graph

    def _set_trap_nodes(self, networkx_graph):
        for node in networkx_graph.nodes():
            networkx_graph.add_node(node, node_type="trap_node", color="b")

    def _remove_horizontal_edges(self, networkx_graph):
        for i in range(0, self.m_extended - self.ion_chain_size_vertical, self.ion_chain_size_vertical):
            for k in range(1, self.ion_chain_size_vertical):
                for j in range(self.n_extended - 1):
                    networkx_graph.remove_edge((i + k, j), (i + k, j + 1))

    def _remove_vertical_edges(self, networkx_graph):
        for i in range(0, self.n_extended - self.ion_chain_size_horizontal, self.ion_chain_size_horizontal):
            for k in range(1, self.ion_chain_size_horizontal):
                for j in range(self.m_extended - 1):
                    networkx_graph.remove_edge((j, i + k), (j + 1, i + k))

    def _remove_horizontal_nodes(self, networkx_graph):
        for i in range(0, self.m_extended - self.ion_chain_size_vertical, self.ion_chain_size_vertical):
            for k in range(1, self.ion_chain_size_vertical):
                for j in range(0, self.n_extended - self.ion_chain_size_horizontal, self.ion_chain_size_horizontal):
                    for s in range(1, self.ion_chain_size_horizontal):
                        networkx_graph.remove_node((i + k, j + s))

    def _set_junction_nodes(self, networkx_graph):
        for i in range(0, self.m_extended, self.ion_chain_size_vertical):
            for j in range(0, self.n_extended, self.ion_chain_size_horizontal):
                networkx_graph.add_node((i, j), node_type="junction_node", color="g")

    def _set_processing_zone(self, networkx_graph):
        # Define the key nodes
        self.exit = (self.m_extended - 1, self.n_extended - 1)
        self.processing_zone = (self.m_extended + self.num_edges - 1, self.n_extended + self.num_edges - 1)
        self.entry = (self.m_extended - 1, 0)

        # differences
        dy_exit = self.exit[1] - self.processing_zone[1]
        dy_entry = self.processing_zone[1] - self.entry[1]

        # Add exit edges
        for i in range(self.num_edges):
            exit_node = (self.exit[0] + (i + 1), self.exit[1] - (i + 1) * dy_exit / self.num_edges)
            networkx_graph.add_node(exit_node, node_type="exit_node", color="y")
            if i == 0:
                previous_exit_node = self.exit

            # changed to exit_node, previous_exit_node
            networkx_graph.add_edge(exit_node, previous_exit_node, edge_type="exit", color="k")

            previous_exit_node = exit_node

        # Add entry edges
        for i in range(self.num_edges):
            entry_node = (self.entry[0] + (i + 1), self.entry[1] + (i + 1) * dy_entry / self.num_edges)
            networkx_graph.add_node(entry_node, node_type="entry_node", color="orange")
            if i == 0:
                previous_entry_node = self.entry

            networkx_graph.add_edge(previous_entry_node, entry_node, edge_type="entry", color="k")

            previous_entry_node = entry_node

        assert exit_node == entry_node, "Exit and entry do not end in same node"
        assert exit_node == self.processing_zone, "Exit and entry do not end in processing zone"

        # Add the processing zone node
        networkx_graph.add_node(self.processing_zone, node_type="processing_zone_node", color="r")

    # plotting
    def plot_state(self, ion_moves, labels, plot_ions=True, show_plot=False):
        # idc_dict = create_idc_dicitonary(nx_G)
        pos = {(x, y): (y, -x) for i, (x, y) in enumerate(list(self.networkx_graph.nodes()))}
        if plot_ions is True:
            pass
            # edge_labels = nx.get_edge_attributes(self.networkx_graph,'ion_chain')
        else:
            edge_labels = {}
            for idc in self.networkx_graph.edges():
                # pass
                edge_labels[idc] = "$e_{%s}$" % get_idx_from_idc(self.idc_dict, idc)

        for edge_idc in self.networkx_graph.edges():
            # color all edges black
            self.networkx_graph.add_edge(edge_idc[0], edge_idc[1], color="k")

            ion_holder = {}
            colors = []
            np.random.seed(0)
            for _ in range(len(ion_moves)):
                r = np.round(np.random.rand(), 1)
                g = np.round(np.random.rand(), 1)
                b = np.round(np.random.rand(), 1)

                colors.append((r, g, b))
            np.random.seed()

            for i, ion_place in enumerate(ion_moves):
                ion_edge_idc = get_idc_from_idx(self.idc_dict, ion_place)
                try:
                    ion_holder[ion_place].append(i)
                except KeyError:
                    ion_holder[ion_place] = [i]
            for i, ion_place in enumerate(ion_moves):
                ion_edge_idc = get_idc_from_idx(self.idc_dict, ion_place)
                self.networkx_graph.add_edge(
                    ion_edge_idc[0], ion_edge_idc[1], ion_chain=ion_holder[ion_place], color=colors[i]
                )

        edge_color = nx.get_edge_attributes(self.networkx_graph, "color").values()
        node_color = list(nx.get_node_attributes(self.networkx_graph, "color").values())
        edge_labels = nx.get_edge_attributes(self.networkx_graph, "ion_chain")

        # plt.figure(figsize=(25, 15))
        plt.figure(figsize=(15, 8))
        nx.draw_networkx(
            self.networkx_graph,
            pos=pos,
            with_labels=True,
            node_size=300,
            node_color=node_color,
            width=8,
            edge_color=edge_color,
            font_size=6,
        )
        nx.draw_networkx_edge_labels(self.networkx_graph, pos, edge_labels)

        # reset edge labels
        for i, ion in enumerate(ion_moves):
            ion_edge_idc = get_idc_from_idx(self.idc_dict, ion)
            self.networkx_graph.add_edge(ion_edge_idc[0], ion_edge_idc[1], ion_chain="", color=colors[i])

        labels0, labels1 = labels
        plt.plot([], [], label=labels0)
        plt.plot([], [], label=labels1)
        plt.legend()

        if show_plot is True:
            plt.show()


GC = GraphCreator(3, 3, 1, 1)
GC.plot_state([0, 1, 2, 3, 4, 5, 6, 7, 8], ("label0", "label1"), show_plot=True)


graph = GC.create_graph()
print(create_idc_dictionary(graph))
