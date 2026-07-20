/*******************************************************************************
 * Copyright (c) 2010, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.graph;

import java.util.Iterator;
import java.util.List;

import org.eclipse.elk.alg.layered.intermediate.LabelAndNodeSizeProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.PortSide;

import java.util.ArrayList;
import java.util.function.Predicate;

/**
 * A port in a layered graph. The position of the port is relative to the upper left corner
 * of the containing node and marks the port's upper left corner. The position where edges are
 * attached is marked by the anchor and is relative to the port's position. The incoming edges
 * and outgoing edges are stored separately, but is is possible to iterate over both types of
 * edges using the concatenation of the corresponding lists.
 * 
 * <p>Port must be used even if the original graph does not reveal them. In this
 * case each edge has dedicated source and target ports, which are used to
 * determine the points where the edge touches the source and target nodes.</p>
 */
public final class LPort extends LShape {

    /** the serial version UID. */
    private static final long serialVersionUID = -3406558719744943360L;
    
    /** a predicate that checks for output ports, that is ports with outgoing edges. */
    public static final Predicate<LPort> OUTPUT_PREDICATE = (port) -> !port.outgoingEdges.isEmpty();
    /** a predicate that checks for input ports, that is ports with incoming edges. */
    public static final Predicate<LPort> INPUT_PREDICATE = (port) -> !port.incomingEdges.isEmpty();
    /** a predicate that checks for north-side ports. */
    public static final Predicate<LPort> NORTH_PREDICATE = (port) -> port.side == PortSide.NORTH;
    /** a predicate that checks for east-side ports. */
    public static final Predicate<LPort> EAST_PREDICATE = (port) -> port.side == PortSide.EAST;
    /** a predicate that checks for south-side ports. */
    public static final Predicate<LPort> SOUTH_PREDICATE = (port) -> port.side == PortSide.SOUTH;
    /** a predicate that checks for west-side ports. */
    public static final Predicate<LPort> WEST_PREDICATE = (port) -> port.side == PortSide.WEST;
    
    /** the owning node. */
    private LNode owner;
    /** the port side. */
    private PortSide side = PortSide.UNDEFINED;
    /** the anchor point position. */
    private final KVector anchor = new KVector();
    /** whether the anchor point position was fixed by the user or not. */
    private boolean explicitlySuppliedPortAnchor = false;
    /** the margin area around this port. */
    private final LMargin margin = new LMargin();
    /** this port's labels. */
    private final List<LLabel> labels = new ArrayList<>(2);

    /** the edges going into the port. */
    private final List<LEdge> incomingEdges = new ArrayList<>(4);
    /** the edges going out of the port. */
    private final List<LEdge> outgoingEdges = new ArrayList<>(4);
    /** All connected edges in a combined iterable. */
    private Iterable<LEdge> connectedEdges = new CombineIter<LEdge>(incomingEdges, outgoingEdges);
    
    /**
     * Whether the original port has edges that connect it to nodes outside of this subgraph. This influences the amount
     * of space we will have to reserve for labels in the {@link LabelAndNodeSizeProcessor}.
     */
    private boolean connectedToExternalNodes = true;
    
    /**
     * Returns the node that owns this port.
     * 
     * @return the owning node
     */
    public LNode getNode() {
        return owner;
    }

    /**
     * Sets the owning node and adds itself to the node's list of ports.
     * If the port was previously in another node, it is removed from that
     * node's list of ports. Be careful not to use this method while
     * iterating through the ports list of the old node nor of the new node,
     * since that could lead to {@link java.util.ConcurrentModificationException}s.
     * 
     * @param node the owner to set
     */
    public void setNode(final LNode node) {
        if (this.owner != null) {
            this.owner.getPorts().remove(this);
        }
        
        this.owner = node;
        
        if (this.owner != null) {
            this.owner.getPorts().add(this);
        }
    }

    /**
     * Returns the node side on which the port is drawn.
     * 
     * @return the side
     */
    public PortSide getSide() {
        return side;
    }

    /**
     * Sets the node side on which the port is drawn.
     * If the port anchor was not already supplied it is set to the default for the corresponding direction.
     * The default position is in the middle of the port in the direction specified by the side to set.
     * 
     * @param theside the side to set
     */
    public void setSide(final PortSide theside) {
        if (theside == null) {
            throw new NullPointerException();
        }
        this.side = theside;
        if (!isExplicitlySuppliedPortAnchor()) {
            switch (this.side) {
            case NORTH:
                // adapt the anchor so edges are attached top center
                getAnchor().x = getSize().x / 2;
                getAnchor().y = 0;
                break;
            case EAST:
                // adapt the anchor so edges are attached center right
                getAnchor().x = getSize().x;
                getAnchor().y = getSize().y / 2;
                break;
            case SOUTH:
                // adapt the anchor so edges are attached bottom center
                getAnchor().x = getSize().x / 2;
                getAnchor().y = getSize().y;
                break;
            case WEST:
                // adapt the anchor so edges are attached center left
                getAnchor().x = 0;
                getAnchor().y = getSize().y / 2;
                break;
            }
        }
    }
    
    /**
     * Returns the anchor position of the port. This is the point where edges should be attached,
     * relative to the port's position. Is automatically changed whenever {@code setSide} is called.
     * 
     * @return the anchor position
     */
    public KVector getAnchor() {
        return anchor;
    }
    
    /**
     * Checks if the anchor position is explicitly set. The algorithm should always check this
     * property before changing the port anchor.
     * 
     * @return {@code true} or {@code false} as the anchor is explicitly supplied, respectively.
     */
    public boolean isExplicitlySuppliedPortAnchor() {
        return explicitlySuppliedPortAnchor;
    }
    
    /**
     * Determines whether the port's anchor is to be considered explicitly supplied and should not be changed.
     * 
     * @param fixed {@code true} if the anchor is explicitly supplied, {@code false} if it is eligible
     *              to be modified.
     */
    public void setExplicitlySuppliedPortAnchor(final boolean fixed) {
        explicitlySuppliedPortAnchor = fixed;
    }
    
    /**
     * Returns the absolute anchor position of the port. This is the point where edges should be
     * attached, relative to the containing graph. This method creates a new vector, so modifying
     * the result will not affect this port in any way.
     * 
     * @return a new vector with the absolute anchor position
     */
    public KVector getAbsoluteAnchor() {
        return KVector.sum(owner.getPosition(), this.getPosition(), anchor);
    }
    
    /**
     * Returns the margin around this port. The margin is typically used to reserve space for the
     * port's labels.
     * 
     * @return the port's margin.
     */
    public LMargin getMargin() {
        return margin;
    }
    
    /**
     * Returns this port's labels.
     * 
     * @return this port's labels.
     */
    public List<LLabel> getLabels() {
        return labels;
    }
    
    /**
     * Returns the name of the port. The name is derived from the text of the first label, if any.
     * 
     * @return the name, or {@code null}
     */
    public String getName() {
        if (!labels.isEmpty()) {
            return labels.get(0).getText();
        }
        return null;
    }
    
    /**
     * Returns this port's degree, that is, the number of edges connected to it.
     * 
     * @return the number of edges connected to this port.
     */
    public int getDegree() {
        return incomingEdges.size() + outgoingEdges.size();
    }
    
    /**
     * Returns the number of incoming edges minus the number of outgoing edges. This
     * is the net flow of the port.
     * 
     * @return the port's net flow.
     */
    public int getNetFlow() {
        return incomingEdges.size() - outgoingEdges.size();
    }

    /**
     * Returns the list of edges going into this port.
     * 
     * @return the incoming edges
     */
    public List<LEdge> getIncomingEdges() {
        return incomingEdges;
    }

    /**
     * Returns the list of edges going out of this port.
     * 
     * @return the outgoing edges
     */
    public List<LEdge> getOutgoingEdges() {
        return outgoingEdges;
    }
    
    /**
     * Returns an iterable over all connected edges, both incoming and outgoing.
     * 
     * @return an iterable over all connected edges.
     */
    public Iterable<LEdge> getConnectedEdges() {
        return connectedEdges;
    }
    
    /**
     * Returns whether the original port has connections to nodes outside of this subgraph.
     */
    public boolean isConnectedToExternalNodes() {
        return connectedToExternalNodes;
    }
    
    /**
     * Returns whether the original port has connections to nodes outside of this subgraph.
     */
    public void setConnectedToExternalNodes(final boolean conn) {
        connectedToExternalNodes = conn;
    }

    /**
     * Returns an iterable over all the port's predecessor ports. Predecessor ports are source
     * ports of incoming edges of this port.
     * 
     * @return an iterable over all predecessor ports.
     */
    public Iterable<LPort> getPredecessorPorts() {
        return new Iterable<LPort>() {
            public Iterator<LPort> iterator() {
                final Iterator<LEdge> edgesIter = incomingEdges.iterator();
                
                return new Iterator<LPort>() {
                    public boolean hasNext() {
                        return edgesIter.hasNext();
                    }
                    public LPort next() {
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
     * Returns an iterable over all the port's successor ports. Successor ports are target
     * ports of outgoing edges of this port.
     * 
     * @return an iterable over all successor ports.
     */
    public Iterable<LPort> getSuccessorPorts() {
        return new Iterable<LPort>() {
            public Iterator<LPort> iterator() {
                final Iterator<LEdge> edgesIter = outgoingEdges.iterator();
                
                return new Iterator<LPort>() {
                    public boolean hasNext() {
                        return edgesIter.hasNext();
                    }
                    public LPort next() {
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
     * Returns an iterable over all connected ports, both predecessors and successors.
     * 
     * @return an iterable over the connected ports
     */
    public Iterable<LPort> getConnectedPorts() {
        return LGraphUtil.concat(getPredecessorPorts(), getSuccessorPorts());
    }
    
    /**
     * Returns the index of the port in the containing node's list of ports. Note
     * that this method has linear running time in the number of ports, so use
     * it with caution.
     * 
     * @return the index of this port, or -1 if the port has no owner
     */
    public int getIndex() {
        if (owner == null) {
            return -1;
        } else {
            return owner.getPorts().indexOf(this);
        }
    }

    @Override
    public String getDesignation() {
        if (!labels.isEmpty() && !(labels.get(0).getText() == null || labels.get(0).getText().isEmpty())) {
            return labels.get(0).getText();
        }
        String id = super.getDesignation();
        if (id != null) {
            return id;
        }
        return Integer.toString(getIndex());
    }

    @Override
    public String toString() {
        StringBuilder result = new StringBuilder();
        result.append("p_").append(getDesignation());
        if (owner != null) {
            result.append("[").append(owner).append("]");
        }
        if (incomingEdges.size() == 1 && outgoingEdges.isEmpty() && incomingEdges.get(0).getSource() != this) {
            LPort source = incomingEdges.get(0).getSource();
            result.append(" << ").append(source.getDesignation());
            result.append("[").append(source.owner).append("]");
        }
        if (incomingEdges.isEmpty() && outgoingEdges.size() == 1 && outgoingEdges.get(0).getTarget() != this) {
            LPort target = outgoingEdges.get(0).getTarget();
            result.append(" >> ").append(target.getDesignation());
            result.append("[").append(target.owner).append("]");
        }
        return result.toString();
    }

    /**
     * Combines two Iterables. Local fast-path replacement for a generic {@code concat} because the
     * lazy implementation here is faster for the small inputs seen in this hot path.
     * 
     * @param <T>
     */
    private static class CombineIter<T> implements Iterable<T> {
        private Iterable<T> firstIterable;
        private Iterable<T> secondIterable;
    
        CombineIter(final Iterable<T> firstIterable, final Iterable<T> secondIterable) {
            this.firstIterable = firstIterable;
            this.secondIterable = secondIterable;
        }
    
        @Override
        public Iterator<T> iterator() {
            return new Iterator<T>() {
                private Iterator<T> firstIterator = firstIterable.iterator();
                private Iterator<T> secondIterator = secondIterable.iterator();
    
                @Override
                public boolean hasNext() {
                    return firstIterator.hasNext() || secondIterator.hasNext();
                }
    
                @Override
                public T next() {
                    if (firstIterator.hasNext()) {
                        return firstIterator.next();
                    } else {
                        return secondIterator.next();
                    }
                }
            };
        }
    }

}
