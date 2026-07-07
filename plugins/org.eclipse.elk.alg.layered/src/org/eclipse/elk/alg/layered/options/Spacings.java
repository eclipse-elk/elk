/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphElement;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;

/**
 * FIXME 
 * Container class for a variety of spacing values that are either specified in the general
 * {@link LayeredOptions} class or ELK Layered's dedicated {@link LayeredOptions} class.
 * 
 * This class allows to either select the recorded spacing values directly or to query for spacing
 * values using one of the convenience methods. The methods do not provide results for every
 * combination of graph elements yet. In case it does not know the answer a
 * {@link UnspecifiedSpacingException} is thrown. In such a case the developer should add the
 * required functionality to this class.
 * 
 * @author uru
 */
public final class Spacings {

    // SUPPRESS CHECKSTYLE NEXT 40 VisibilityModifier
    
    private LGraph graph;
    
    private  IProperty<Double>[][] nodeTypeSpacingOptionsHorizontal;
    private  IProperty<Double>[][] nodeTypeSpacingOptionsVertical;
    
    
    /**
     * @param graph
     *            the {@link LGraph} for which to record the spacing values.
     */
    @SuppressWarnings("unchecked")
    public Spacings(final LGraph graph) {
        this.graph = graph;

        // pre calculate the spacings between pairs of node types
        int n = NodeType.values().length;
        nodeTypeSpacingOptionsHorizontal = new IProperty[n][n];
        nodeTypeSpacingOptionsVertical = new IProperty[n][n];
        precalculateNodeTypeSpacings();
    }
    
    private void precalculateNodeTypeSpacings() {
        
        // normal
        nodeTypeSpacing(NodeType.NORMAL, 
                LayeredOptions.SPACING_NODE_NODE,
                LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS);
        nodeTypeSpacing(NodeType.NORMAL, NodeType.LONG_EDGE,
                LayeredOptions.SPACING_EDGE_NODE,
                LayeredOptions.SPACING_EDGE_NODE_BETWEEN_LAYERS);
        nodeTypeSpacing(NodeType.NORMAL, NodeType.NORTH_SOUTH_PORT,
                LayeredOptions.SPACING_EDGE_NODE);
        nodeTypeSpacing(NodeType.NORMAL, NodeType.EXTERNAL_PORT,
                LayeredOptions.SPACING_EDGE_NODE); // TODO
        nodeTypeSpacing(NodeType.NORMAL, NodeType.LABEL,
                LayeredOptions.SPACING_NODE_NODE,
                LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS);
        
        // longedge
        nodeTypeSpacing(NodeType.LONG_EDGE, 
                LayeredOptions.SPACING_EDGE_EDGE,
                LayeredOptions.SPACING_EDGE_EDGE_BETWEEN_LAYERS);
        nodeTypeSpacing(NodeType.LONG_EDGE, NodeType.NORTH_SOUTH_PORT, 
                LayeredOptions.SPACING_EDGE_EDGE);
        nodeTypeSpacing(NodeType.LONG_EDGE, NodeType.EXTERNAL_PORT, 
                LayeredOptions.SPACING_EDGE_EDGE); // TODO
        nodeTypeSpacing(NodeType.LONG_EDGE, NodeType.LABEL, 
                LayeredOptions.SPACING_EDGE_NODE,
                LayeredOptions.SPACING_EDGE_NODE_BETWEEN_LAYERS);

        // northsouth
        nodeTypeSpacing(NodeType.NORTH_SOUTH_PORT, 
                LayeredOptions.SPACING_EDGE_EDGE);
        nodeTypeSpacing(NodeType.NORTH_SOUTH_PORT, NodeType.EXTERNAL_PORT, 
                LayeredOptions.SPACING_EDGE_EDGE); // TODO
        nodeTypeSpacing(NodeType.NORTH_SOUTH_PORT, NodeType.LABEL, 
                LayeredOptions.SPACING_LABEL_NODE);

        // external
        nodeTypeSpacing(NodeType.EXTERNAL_PORT, 
                LayeredOptions.SPACING_PORT_PORT);
        nodeTypeSpacing(NodeType.EXTERNAL_PORT, NodeType.LABEL, 
                LayeredOptions.SPACING_LABEL_PORT_VERTICAL, LayeredOptions.SPACING_LABEL_PORT_HORIZONTAL);
        
        // label
        nodeTypeSpacing(NodeType.LABEL,
                LayeredOptions.SPACING_EDGE_EDGE,
                LayeredOptions.SPACING_EDGE_EDGE);
        
        // breaking points
        nodeTypeSpacing(NodeType.BREAKING_POINT, 
                LayeredOptions.SPACING_EDGE_EDGE,
                LayeredOptions.SPACING_EDGE_EDGE_BETWEEN_LAYERS);
        nodeTypeSpacing(NodeType.BREAKING_POINT, NodeType.NORMAL, 
                LayeredOptions.SPACING_EDGE_NODE,
                LayeredOptions.SPACING_EDGE_NODE_BETWEEN_LAYERS);
        nodeTypeSpacing(NodeType.BREAKING_POINT, NodeType.LABEL,  
                LayeredOptions.SPACING_EDGE_NODE, 
                LayeredOptions.SPACING_EDGE_NODE_BETWEEN_LAYERS); 
        nodeTypeSpacing(NodeType.BREAKING_POINT, NodeType.LONG_EDGE, 
                LayeredOptions.SPACING_EDGE_NODE,
                LayeredOptions.SPACING_EDGE_NODE_BETWEEN_LAYERS);
        
    }

    private void nodeTypeSpacing(final NodeType nt, final IProperty<Double> spacing) {
        nodeTypeSpacingOptionsVertical[nt.ordinal()][nt.ordinal()] = spacing;
    }

    private void nodeTypeSpacing(final NodeType nt, final IProperty<Double> spacingVert,
            final IProperty<Double> spacingHorz) {
        nodeTypeSpacingOptionsVertical[nt.ordinal()][nt.ordinal()] = spacingVert;

        nodeTypeSpacingOptionsHorizontal[nt.ordinal()][nt.ordinal()] = spacingHorz;
    }

    private void nodeTypeSpacing(final NodeType n1, final NodeType n2, final IProperty<Double> spacing) {
        nodeTypeSpacingOptionsVertical[n1.ordinal()][n2.ordinal()] = spacing;
        nodeTypeSpacingOptionsVertical[n2.ordinal()][n1.ordinal()] = spacing;
    }

    private void nodeTypeSpacing(final NodeType n1, final NodeType n2, final IProperty<Double> spacingVert,
            final IProperty<Double> spacingHorz) {
        nodeTypeSpacingOptionsVertical[n1.ordinal()][n2.ordinal()] = spacingVert;
        nodeTypeSpacingOptionsVertical[n2.ordinal()][n1.ordinal()] = spacingVert;

        nodeTypeSpacingOptionsHorizontal[n1.ordinal()][n2.ordinal()] = spacingHorz;
        nodeTypeSpacingOptionsHorizontal[n2.ordinal()][n1.ordinal()] = spacingHorz;
    }

    /**
     * @param e1
     *            a graph element
     * @param e2
     *            another graph element
     * @return the horizontal spacing to be preserved between {@code e1} and {@code e2}
     */
    public double getHorizontalSpacing(final LGraphElement e1, final LGraphElement e2) {
        if (e1 instanceof LNode && e2 instanceof LNode) {
            return getHorizontalSpacing((LNode) e1, (LNode) e2);
        }

        throw new UnspecifiedSpacingException();
    }

    /**
     * @param n1
     *            a node
     * @param n2
     *            another node
     * @return the spacing to be preserved between {@code n1} and {@code n2}
     */
    public double getHorizontalSpacing(final LNode n1, final LNode n2) {
        return getLocalSpacing(n1, n2, nodeTypeSpacingOptionsHorizontal);
    }
    
    /**
     * @param nt1
     *            a node type
     * @param nt2
     *            another node type
     * @return the spacing to be preserved between {@code nt1} and {@code nt2}
     */
    public double getHorizontalSpacing(final NodeType nt1, final NodeType nt2) {
        return getLocalSpacing(nt1, nt2, nodeTypeSpacingOptionsHorizontal);
    }
    
    /**
     * @param n1
     *            a node
     * @param n2
     *            another node
     * @return the spacing to be preserved between {@code n1} and {@code n2}
     */
    public double getVerticalSpacing(final LNode n1, final LNode n2) {
        return getLocalSpacing(n1, n2, nodeTypeSpacingOptionsVertical);
    }

    /**
     * @param nt1
     *            a node
     * @param nt2
     *            another node
     * @return the spacing to be preserved between {@code n1} and {@code n2}
     */
    public double getVerticalSpacing(final NodeType nt1, final NodeType nt2) {
        return getLocalSpacing(nt1, nt2, nodeTypeSpacingOptionsVertical);
    }

    private double getLocalSpacing(final LNode n1, final LNode n2, final IProperty<Double>[][] nodeTypeSpacingMapping) {
        NodeType t1 = n1.getType();
        NodeType t2 = n2.getType();
        IProperty<Double> layoutOption = nodeTypeSpacingMapping[t1.ordinal()][t2.ordinal()];

        // get the spacing value for the first node
        Double s1 = getIndividualOrDefault(n1, layoutOption); 
        Double s2 = getIndividualOrDefault(n2, layoutOption);
        
        return Math.max(s1, s2);
    }

    private double getLocalSpacing(final NodeType nt1, final NodeType nt2,
            final IProperty<Double>[][] nodeTypeSpacingMapping) {
        IProperty<Double> layoutOption = nodeTypeSpacingMapping[nt1.ordinal()][nt2.ordinal()];
        return graph.getProperty(layoutOption);
    }

    /**
     * Returns the value of the given property as it applies to the given node. First checks whether an individual
     * override is set on the node that has the given property configured. If so, the configured value is returned.
     * Otherwise, the node's parent node, if any, is queried.
     * 
     * @param node
     *            the node whose property value to return.
     * @param property
     *            the property.
     * @return the individual override for the property or the default value inherited by the parent node.
     */
    public static <T> T getIndividualOrDefault(final LNode node, final IProperty<T> property) {
        T result = null;
        // check for individual value
        if (node.hasProperty(LayeredOptions.SPACING_INDIVIDUAL)) {
            IPropertyHolder individualSpacings = node.getProperty(LayeredOptions.SPACING_INDIVIDUAL);
            if (individualSpacings.hasProperty(property)) {
                result = individualSpacings.getProperty(property);
            }
        }
        // use the common value
        if (result == null) {
            result = node.getGraph().getProperty(property);
        }
        
        return result;
    }

    /**
     * Dedicated exception indicating that no spacing value could be determined for a certain set of
     * graph elements. This is probably due to a programming error.
     */
    public static class UnspecifiedSpacingException extends RuntimeException {

        private static final long serialVersionUID = 1609767701465615319L;

        // SUPPRESS CHECKSTYLE NEXT 4 Javadoc
        public UnspecifiedSpacingException() {
        }

        public UnspecifiedSpacingException(final String msg) {
            super(msg);
        }
    }
}
