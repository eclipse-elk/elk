/*******************************************************************************
 * Copyright (c) 2014, 2019 Kiel University and others.
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

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.Predicate;
import java.util.ArrayList;

import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopEdge;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopPort;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.LabelSide;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.adapters.GraphAdapters.EdgeAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphElementAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.LabelAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.NodeAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.PortAdapter;
import org.eclipse.elk.graph.properties.IProperty;


/**
 * Provides implementations of the {@link org.eclipse.elk.core.util.adapters.GraphAdapters
 * GraphAdapters} interfaces for the LGraph. The adapted graph can then be fed to ELK's node size
 * calculation code, for example. To obtain an adapter for an {@link LGraph}, simply call
 * {@link #adapt(LGraph)}.
 */
public final class LGraphAdapters {

    private LGraphAdapters() {
        throw new IllegalStateException("Private constructor, not to be instantiated!");
    }
    
    
    /**
     * Adapts the given {@link LGraph}.
     * 
     * @param graph
     *            the graph that should be wrapped in an adapter.
     * @return an {@link LGraphAdapter} for the passed graph.
     */
    public static LGraphAdapter adapt(final LGraph graph) {
        return adapt(graph, false);
    }
    
    /**
     * Adapts the given {@link LGraph}. Transparently provides access to edges connected to north/south port dummies.
     *
     * @param graph
     *            the graph that should be wrapped in an adapter.
     * @param transparentNorthSouthEdges
     *            {@code true} if edges connected to north south port dummies should appear to be
     *            directly connected to their original north/south ports. This effectively makes the
     *            north/south port dummies "transparent" in the sense that edges connected to them
     *            appear to be connected to the original connection points.
     * @return an {@link LGraphAdapter} for the passed graph.
     */
    public static LGraphAdapter adapt(final LGraph graph, final boolean transparentNorthSouthEdges) {
        return new LGraphAdapter(graph, transparentNorthSouthEdges, false, n -> true);
     }
    
    /**
     * Adapts the given {@link LGraph}. Transparently provides access to edges connected to
     * north/south port dummies and allows to filter out nodes that client code shouldn't see.
     *
     * @param graph
     *            the graph that should be wrapped in an adapter.
     * @param transparentNorthSouthEdges
     *            {@code true} if edges connected to north south port dummies should appear to be
     *            directly connected to their original north/south ports. This effectively makes the
     *            north/south port dummies "transparent" in the sense that edges connected to them
     *            appear to be connected to the original connection points.
     * @param transparentCommentNodes
     *            {@code true} if attached comment nodes temporarily removed from the graph should
     *            be included in the adapted graph as if they were included.
     * @param nodeFilter
     *            predicate that returns {@code true} for a node that should be visible in the adapted
     *            graph and {@code false} for one that should not. Usual application is to hide
     *            certain kinds of dummy nodes.
     * @return an {@link LGraphAdapter} for the passed graph.
     */
    public static LGraphAdapter adapt(final LGraph graph, final boolean transparentNorthSouthEdges,
            final boolean transparentCommentNodes, final Predicate<LNode> nodeFilter) {
        
        return new LGraphAdapter(graph, transparentNorthSouthEdges, transparentCommentNodes, nodeFilter);
    }
    
    /**
     * Adapts the given single node. Transparently provides access to edges connected to north/south port dummies.
     *
     * @param node
     *            the node that should be wrapped in an adapter.
     * @param transparentNorthSouthEdges
     *            {@code true} if edges connected to north south port dummies should appear to be
     *            directly connected to their original north/south ports. This effectively makes the
     *            north/south port dummies "transparent" in the sense that edges connected to them
     *            appear to be connected to the original connection points.
     * @return an {@link LNodeAdapter} for the passed node.
     */
    public static LNodeAdapter adapt(final LNode node, final boolean transparentNorthSouthEdges) {
        return new LNodeAdapter(null, node, transparentNorthSouthEdges);
    }
    
    /**
     * Adapts the given {@link LLabel}.
     * 
     * @param label
     *            the label to wrap an adapter around.
     * @return an {@link LLabelAdapter} for the passed label.
     */
    public static LLabelAdapter adapt(final LLabel label) {
        return new LLabelAdapter(label);
    }
    
    
    /**
     * Basic base class for adapters that adapt {@link LGraphElement}s.
     * @param <T> type of the adapted element.
     */
    private abstract static class AbstractLShapeAdapter<T extends LShape> implements GraphElementAdapter<T> {

        // CHECKSTYLEOFF VisibilityModifier
        /** The wrapped element. */
        protected T element;
        // CHECKSTYLEON VisibilityModifier
        
        
        /**
         * Creates a new adapter for the given element.
         * 
         * @param element
         *            the element to be adapted.
         */
        AbstractLShapeAdapter(final T element) {
            this.element = element;
        }
        

        @Override
        public KVector getSize() {
            return element.getSize();
        }

        @Override
        public void setSize(final KVector size) {
            element.getSize().x = size.x;
            element.getSize().y = size.y;
        }

        @Override
        public KVector getPosition() {
            return element.getPosition();
        }

        @Override
        public void setPosition(final KVector pos) {
            element.getPosition().x = pos.x;
            element.getPosition().y = pos.y;
        }

        @Override
        public <P> P getProperty(final IProperty<P> prop) {
            return element.getProperty(prop);
        }

        @Override
        public <P> boolean hasProperty(final IProperty<P> prop) {
            return element.hasProperty(prop);
        }

        @Override
        public int getVolatileId() {
            return element.id;
        }

        @Override
        public void setVolatileId(final int volatileId) {
            element.id = volatileId;
        }
    }

    /**
     * Adapter for {@link LGraph}s.
     */
    public static final class LGraphAdapter implements GraphAdapter<LGraph> {
        // CHECKSTYLEOFF VisibilityModifier
        /** The wrapped element. */
        protected final LGraph element;
        // CHECKSTYLEON VisibilityModifier
        /** List of cached node adapters. */
        private List<NodeAdapter<?>> nodeAdapters = null;
        /**
         * Whether to simulate that edges are directly connected to north south ports instead of to
         * north/south port dummies.
         */
        private final boolean transparentNorthSouthEdges;
        /**
         * Whether to include attached comment nodes currently not part of the graph.
         */
        private final boolean transparentCommentNodes;
        /** Predicate that decides which of our children we make visible to clients. */
        private final Predicate<LNode> nodeFilter;

        /**
         * Creates a new adapter for the given graph.
         * 
         * @param element
         *            the graph to be adapted.
         * @param transparentNorthSouthEdges
         *            whether to simulate that edges are directly connected to north south ports 
         *            instead of to north/south port dummies.
         */
        private LGraphAdapter(final LGraph element, final boolean transparentNorthSouthEdges,
                final boolean transparentCommentNodes, final Predicate<LNode> nodeFilter) {
            
            this.element = element;
            this.transparentNorthSouthEdges = transparentNorthSouthEdges;
            this.transparentCommentNodes = transparentCommentNodes;
            this.nodeFilter = nodeFilter;
        }

        @Override
        public KVector getSize() {
            return element.getSize();
        }

        @Override
        public void setSize(final KVector size) {
            element.getSize().x = size.x;
            element.getSize().y = size.y;
        }

        @Override
        public KVector getPosition() {
            throw new UnsupportedOperationException("Not supported by LGraph");
        }

        @Override
        public void setPosition(final KVector pos) {
            throw new UnsupportedOperationException("Not supported by LGraph");
        }

        @Override
        public <P> P getProperty(final IProperty<P> prop) {
            return element.getProperty(prop);
        }

        @Override
        public <P> boolean hasProperty(final IProperty<P> prop) {
            return element.hasProperty(prop);
        }

        @Override
        public Iterable<NodeAdapter<?>> getNodes() {
            if (nodeAdapters == null) {
                nodeAdapters = new ArrayList<>();
                // We completely ignore layerless nodes here since they are currently not of interest
                // to anyone using these adapters
                for (Layer l : element.getLayers()) {
                    for (LNode n : l.getNodes()) {
                        if (nodeFilter.test(n)) {
                            nodeAdapters.add(new LNodeAdapter(this, n, transparentNorthSouthEdges));
                            
                            if (transparentCommentNodes) {
                                if (n.hasProperty(InternalProperties.TOP_COMMENTS)) {
                                    for (LNode comment : n.getProperty(InternalProperties.TOP_COMMENTS)) {
                                        nodeAdapters.add(new LNodeAdapter(this, comment, false));
                                    }
                                }
                                
                                if (n.hasProperty(InternalProperties.BOTTOM_COMMENTS)) {
                                    for (LNode comment : n.getProperty(InternalProperties.BOTTOM_COMMENTS)) {
                                        nodeAdapters.add(new LNodeAdapter(this, comment, false));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return nodeAdapters;
        }

        @Override
        public int getVolatileId() {
            return element.id;
        }

        @Override
        public void setVolatileId(final int volatileId) {
            element.id = volatileId;
        }
    }

    /**
     * Adapter for {@link LNode}s.
     */
    static final class LNodeAdapter extends AbstractLShapeAdapter<LNode> implements NodeAdapter<LNode> {
        /** The parent graph adapter that created this thing. */
        private LGraphAdapter parentGraphAdapter = null;
        /** List of cached label adapters. */
        private List<LabelAdapter<?>> labelAdapters = null;
        /** List of cached port adapters. */
        private List<PortAdapter<?>> portAdapters = null;
        /**
         * Whether to simulate that edges are directly connected to north south ports instead of to
         * north/south port dummies.
         */
        private boolean transparentNorthSouthEdges;

        /**
         * Creates a new adapter for the given node.
         * 
         * @param parent
         *            the graph adapter that created us.
         * @param element
         *            the node to be adapted.
         * @param transparentNorthSouthEdges
         *            whether to simulate that edges are directly connected to north south ports
         *            instead of to north/south port dummies.
         */
        LNodeAdapter(final LGraphAdapter parent, final LNode element, final boolean transparentNorthSouthEdges) {
            super(element);
            this.parentGraphAdapter = parent;
            this.transparentNorthSouthEdges = transparentNorthSouthEdges;
        }
        

        @Override
        public GraphAdapter<?> getGraph() {
            return parentGraphAdapter;
        }

        @Override
        public Iterable<LabelAdapter<?>> getLabels() {
            if (labelAdapters == null) {
                labelAdapters = new ArrayList<>(element.getLabels().size());
                for (LLabel l : element.getLabels()) {
                    labelAdapters.add(new LLabelAdapter(l));
                }
            }
            return labelAdapters;
        }

        @Override
        public Iterable<PortAdapter<?>> getPorts() {
            if (portAdapters == null) {
                portAdapters = new ArrayList<>(element.getPorts().size());
                for (LPort p : element.getPorts()) {
                    portAdapters.add(new LPortAdapter(p, transparentNorthSouthEdges));
                }
            }
            return portAdapters;
        }

        @Override
        public Iterable<EdgeAdapter<?>> getIncomingEdges() {
            // we have no directly connected edges
            return Collections.emptyList();
        }

        @Override
        public Iterable<EdgeAdapter<?>> getOutgoingEdges() {
         // we have no directly connected edges
            return Collections.emptyList();
        }

        @Override
        public void sortPortList() {
            sortPortList(DEFAULT_PORTLIST_SORTER);
        }

        @Override
        @SuppressWarnings("unchecked")
        public void sortPortList(final Comparator<?> comparator) {
            if (element.getProperty(LayeredOptions.PORT_CONSTRAINTS).isOrderFixed()) {
                // We need to sort the port list accordingly
                Collections.sort(element.getPorts(), (Comparator<LPort>) comparator);
            }
        }

        @Override
        public boolean isCompoundNode() {
            return element.getProperty(InternalProperties.COMPOUND_NODE);
        }

        @Override
        public ElkPadding getPadding() {
            LPadding lPadding = element.getPadding();
            return new ElkPadding(lPadding.top, lPadding.right, lPadding.bottom, lPadding.left);
        }

        @Override
        public void setPadding(final ElkPadding padding) {
            element.getPadding().left = padding.left;
            element.getPadding().top = padding.top;
            element.getPadding().right = padding.right;
            element.getPadding().bottom = padding.bottom;
        }

        @Override
        public ElkMargin getMargin() {
            LMargin lmargins = element.getMargin();
            return new ElkMargin(lmargins.top, lmargins.right, lmargins.bottom, lmargins.left);
        }

        @Override
        public void setMargin(final ElkMargin margin) {
            element.getMargin().left = margin.left;
            element.getMargin().top = margin.top;
            element.getMargin().right = margin.right;
            element.getMargin().bottom = margin.bottom;
        }
    }

    /**
     * Adapter for {@link LPort}s.
     */
    static final class LPortAdapter extends AbstractLShapeAdapter<LPort> implements PortAdapter<LPort> {
        /** List of cached label adapters. */
        private List<LabelAdapter<?>> labelAdapters = null;
        /** List of cached edge adapters for incoming edges. */
        private List<EdgeAdapter<?>> incomingEdgeAdapters = null;
        /** List of cached edge adapters for outgoing edges. */
        private List<EdgeAdapter<?>> outgoingEdgeAdapters = null;
        /**
         * Whether to simulate that edges are directly connected to north south ports instead of to
         * north/south port dummies.
         */
        private boolean transparentNorthSouthEdges;

        /**
         * Creates a new adapter for the given port.
         * 
         * @param element
         *            the port to be adapted.
         * @param transparentNorthSouthEdges
         *            whether to simulate that edges are directly connected to north south ports
         *            instead of to north/south port dummies.
         */
        LPortAdapter(final LPort element, final boolean transparentNorthSouthEdges) {
            super(element);
            this.transparentNorthSouthEdges = transparentNorthSouthEdges;
        }
        

        @Override
        public PortSide getSide() {
            return element.getSide();
        }

        @Override
        public Iterable<LabelAdapter<?>> getLabels() {
            if (labelAdapters == null) {
                labelAdapters = new ArrayList<>(element.getLabels().size());
                for (LLabel l : element.getLabels()) {
                    labelAdapters.add(new LLabelAdapter(l));
                }
            }
            return labelAdapters;
        }

        @Override
        public ElkMargin getMargin() {
            LMargin lmargins = element.getMargin();
            return new ElkMargin(lmargins.top, lmargins.right, lmargins.bottom, lmargins.left);
        }

        @Override
        public void setMargin(final ElkMargin margin) {
            element.getMargin().left = margin.left;
            element.getMargin().top = margin.top;
            element.getMargin().right = margin.right;
            element.getMargin().bottom = margin.bottom;
        }

        @Override
        public Iterable<EdgeAdapter<?>> getIncomingEdges() {
            if (transparentNorthSouthEdges && element.getNode().getType() == NodeType.NORTH_SOUTH_PORT) {
                return Collections.emptyList();
            } else if (incomingEdgeAdapters == null) {
                incomingEdgeAdapters = new ArrayList<>();
                for (LEdge e : element.getIncomingEdges()) {
                    incomingEdgeAdapters.add(new LEdgeAdapter(e));
                }
                if (transparentNorthSouthEdges) {
                    final LNode portDummy = element.getProperty(InternalProperties.PORT_DUMMY);
                    if (portDummy != null) {
                        for (LEdge e : portDummy.getIncomingEdges()) {
                            incomingEdgeAdapters.add(new LEdgeAdapter(e));
                       }
                    }
                }
                // Add the incoming edge from the self loop holder if asked for one.
                if (element.getNode().hasProperty(InternalProperties.SELF_LOOP_HOLDER)) {
                    SelfLoopHolder slh = element.getNode().getProperty(InternalProperties.SELF_LOOP_HOLDER);
                    SelfLoopPort slp = slh.getSLPortMap().get(element);
                    if (slp != null) {
                        for (SelfLoopEdge sle : slp.getIncomingSLEdges()) {
                            incomingEdgeAdapters.add(new LEdgeAdapter(sle.getLEdge()));
                        }                        
                    }
                }
            }
            return incomingEdgeAdapters;
        }

        @Override
        public Iterable<EdgeAdapter<?>> getOutgoingEdges() {
            if (transparentNorthSouthEdges && element.getNode().getType() == NodeType.NORTH_SOUTH_PORT) {
                return Collections.emptyList();
            } else if (outgoingEdgeAdapters == null) {
                outgoingEdgeAdapters = new ArrayList<>();
                for (LEdge e : element.getOutgoingEdges()) {
                    outgoingEdgeAdapters.add(new LEdgeAdapter(e));
                }
                if (transparentNorthSouthEdges) {
                    final LNode portDummy = element.getProperty(InternalProperties.PORT_DUMMY);
                    if (portDummy != null) {
                        for (LEdge e : portDummy.getOutgoingEdges()) {
                            outgoingEdgeAdapters.add(new LEdgeAdapter(e));
                        }
                    }
                }
                // Add the outgoing edge from the self loop holder if asked for one.
                if (element.getNode().hasProperty(InternalProperties.SELF_LOOP_HOLDER)) {
                    SelfLoopHolder slh = element.getNode().getProperty(InternalProperties.SELF_LOOP_HOLDER);
                    SelfLoopPort slp = slh.getSLPortMap().get(element);
                    if (slp != null) {
                        for (SelfLoopEdge sle : slp.getOutgoingSLEdges()) {
                            outgoingEdgeAdapters.add(new LEdgeAdapter(sle.getLEdge()));
                        }                        
                    }
                }
            }
            return outgoingEdgeAdapters;
        }

        @Override
        public boolean hasCompoundConnections() {
            return element.getProperty(InternalProperties.INSIDE_CONNECTIONS);
        }
    }

    /**
     * Adapter for {@link LLabel}s.
     */
    static final class LLabelAdapter extends AbstractLShapeAdapter<LLabel> implements
        LabelAdapter<LLabel> {

        /**
         * Creates a new adapter for the given label.
         * 
         * @param element
         *            the label to be adapted.
         */
        LLabelAdapter(final LLabel element) {
            super(element);
        }

        @Override
        public LabelSide getSide() {
            return element.getProperty(LabelSide.LABEL_SIDE);
        }
        
        @Override
        public String getText() {
            return element.getText();
        }
        
    }

    /**
     * Adapter for {@link LEdge}s.
     */
    static final class LEdgeAdapter implements EdgeAdapter<LEdge> {
        /** The wrapped edge. */
        private LEdge element;
        /** List of cached label adapters. */
        private List<LabelAdapter<?>> labelAdapters = null;

        
        /**
         * Creates a new adapter for the given edge.
         * 
         * @param edge
         *            the edge to adapt.
         */
        LEdgeAdapter(final LEdge edge) {
            this.element = edge;
        }
        

        @Override
        public Iterable<LabelAdapter<?>> getLabels() {
            if (labelAdapters == null) {
                labelAdapters = new ArrayList<>(element.getLabels().size());
                for (LLabel l : element.getLabels()) {
                    labelAdapters.add(new LLabelAdapter(l));
                }
            }
            return labelAdapters;
        }
    }

    /**
     * A comparer for ports. Ports are sorted by side (north, east, south, west) in
     * clockwise order, beginning at the top left corner.
     */
    public static final PortComparator DEFAULT_PORTLIST_SORTER = new PortComparator();
    
    /**
     * A comparer for ports. Ports are sorted by side (north, east, south, west) in
     * clockwise order, beginning at the top left corner.
     */
    public static class PortComparator implements Comparator<LPort> {
        @Override
        public int compare(final LPort port1, final LPort port2) {
            int ordinalDifference = port1.getSide().ordinal() - port2.getSide().ordinal();
            
            // Sort by side first
            if (ordinalDifference != 0) {
                return ordinalDifference;
            }
            
            // In case of equal sides, sort by port index property
            Integer index1 = port1.getProperty(LayeredOptions.PORT_INDEX);
            Integer index2 = port2.getProperty(LayeredOptions.PORT_INDEX);
            if (index1 != null && index2 != null) {
                int indexDifference = index1 - index2;
                if (indexDifference != 0) {
                    return indexDifference;
                }
            }
            
            // In case of equal index, sort by position
            switch (port1.getSide()) {
            case NORTH:
                // Compare x coordinates
                return Double.compare(port1.getPosition().x, port2.getPosition().x);
            
            case EAST:
                // Compare y coordinates
                return Double.compare(port1.getPosition().y, port2.getPosition().y);
            
            case SOUTH:
                // Compare x coordinates in reversed order
                return Double.compare(port2.getPosition().x, port1.getPosition().x);
            
            case WEST:
                // Compare y coordinates in reversed order
                return Double.compare(port2.getPosition().y, port1.getPosition().y);
                
            default:
                // Port sides should not be undefined
                throw new IllegalStateException("Port side is undefined");
            }
        }
    }
}
