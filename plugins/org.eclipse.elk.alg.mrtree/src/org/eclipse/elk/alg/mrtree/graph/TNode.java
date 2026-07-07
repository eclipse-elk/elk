/*******************************************************************************
 * Copyright (c) 2013, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.graph;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.alg.mrtree.options.InternalProperties;

/**
 * A node in the T graph. Some properties are maybe null.
 * 
 * @author sor
 * @author sgu
 */
public class TNode extends TShape {

    /** the serial version UID. */
    private static final long serialVersionUID = 1L;

    /** the graph it belongs to. */
    private TGraph graph;

    /** the node label. */
    private String label;

    /** List of outgoing edges. */
    private LinkedList<TEdge> outgoingEdges = new LinkedList<TEdge>();

    /** List of incoming edges. */
    private LinkedList<TEdge> incomingEdges = new LinkedList<TEdge>();

    // CONSTRUCTORS

    /**
     * Create a new node with given label.
     * 
     * @param id
     *            a unique identification
     * @param graph
     *            the graph it belong s to
     * @param label
     *            the label
     */
    public TNode(final int id, final TGraph graph, final String label) {
        this(id, graph);
        this.label = label;
    }

    /**
     * Create a new node a label.
     * 
     * @param id
     *            a unique identification
     * @param graph
     *            the graph it belong s to
     */
    public TNode(final int id, final TGraph graph) {
        super(id);
        this.graph = graph;
    }

    // GETTERS

    @Override
    public String toString() {
        if (label == null || label.length() == 0) {
            return "n_" + id;
        } else {
            return "n_" + label;
        }
    }

    /**
     * Returns the label text of this node.
     * 
     * @return the label
     */
    public String getLabel() {
        return label;
    }

    /**
     * Returns an iterable over all the parents. Parents are sources of incoming edges of this node.
     * 
     * @return an iterable over all parents.
     */
    public Iterable<TNode> getParents() {
        return new Iterable<TNode>() {
            public Iterator<TNode> iterator() {
                final Iterator<TEdge> edgesIter = getIncomingEdges().iterator();

                return new Iterator<TNode>() {
                    public boolean hasNext() {
                        return edgesIter.hasNext();
                    }

                    public TNode next() {
                        return edgesIter.next().getSource();
                    }

                    public void remove() {
                        edgesIter.remove();
                    }
                };
            }

        };
    }

    /**
     * Returns the first parent, by taking it from outgoing edges.
     * 
     * @return a single parent or null if there is none
     */
    public TNode getParent() {
        List<TEdge> edges = getIncomingEdges();
        if (edges.isEmpty()) {
            return null;
        }
        return edges.get(0).getSource();
    }

    /**
     * Returns the list of children, creating it from outgoing edges. Children are targets of
     * outgoing edges of this node.
     * 
     * @return the list of children
     */
    public LinkedList<TNode> getChildrenCopy() {
        LinkedList<TNode> children = new LinkedList<TNode>();
        for (TEdge iEdge : getOutgoingEdges()) {
            children.add(iEdge.getTarget());
        }
        return children;
    }

    /**
     * Returns an iterable over all the children. Children are targets of outgoing edges of this
     * node.
     * 
     * @return an iterable over all children.
     */
    public Iterable<TNode> getChildren() {
        return new Iterable<TNode>() {
            public Iterator<TNode> iterator() {
                final Iterator<TEdge> edgesIter = getOutgoingEdges().iterator();

                return new Iterator<TNode>() {
                    public boolean hasNext() {
                        return edgesIter.hasNext();
                    }

                    public TNode next() {
                        return edgesIter.next().getTarget();
                    }

                    public void remove() {
                        edgesIter.remove();
                    }
                };
            }

        };
    }

    /**
     * Returns the list of outgoing edges, creating it if necessary.
     * 
     * @return list of outgoing edges
     */
    public List<TEdge> getOutgoingEdges() {
        if (outgoingEdges == null) {
            outgoingEdges = new LinkedList<TEdge>();
        }
        return outgoingEdges;
    }

    /**
     * Returns the list of incoming edges, creating it if necessary.
     * 
     * @return list of incoming edge
     */
    public List<TEdge> getIncomingEdges() {
        if (incomingEdges == null) {
            incomingEdges = new LinkedList<TEdge>();
        }
        return incomingEdges;
    }

    // SETTERS

    /**
     * Set the label text of this node.
     * 
     * @param label the new label text
     */
    public void setLabel(final String label) {
        this.label = label;
    }

    /**
     * Add a node to the list of children, by adding a appropriate dummy edge to the graph.
     * 
     * @param child
     *            the child to be added
     */
    public void addChild(final TNode child) {
        TEdge newEdge = new TEdge(this, child);
        newEdge.setProperty(InternalProperties.DUMMY, true);
        graph.getEdges().add(newEdge);
        getOutgoingEdges().add(newEdge);
        child.getIncomingEdges().add(newEdge);
    }

    // ATTRIBUTES

    /**
     * Returns whether this node is a leaf.
     * 
     * @return true if node is a leaf
     */
    public boolean isLeaf() {
        return getOutgoingEdges().isEmpty();
    }
}
