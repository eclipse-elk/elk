/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPadding;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.ContentAlignment;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.options.SizeOptions;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * This processor resizes a child graph to fit the parent node. It must be run as the last non-hierarchical processor in
 * a hierarchical graph.
 * 
 * <dl>
 * <dt>Precondition:</dt>
 * <dd>graph with layout completed</dd>
 * <dt>Postcondition:</dt>
 * <dd>Graph is resized to fit parent node</dd>
 * <dt>Slots:</dt>
 * <dd>After phase 5.</dd>
 */
public class HierarchicalNodeResizingProcessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Resize child graph to fit parent.", 1);
        for (Layer layer : graph) {
            graph.getLayerlessNodes().addAll(layer.getNodes());
            layer.getNodes().clear();
        }
        for (LNode node : graph.getLayerlessNodes()) {
            node.setLayer(null);
        }
        graph.getLayers().clear();
        resizeGraph(graph);
        if (isNested(graph)) {
            graphLayoutToNode(graph.getParentNode(), graph);
        }
        progressMonitor.done();
    }

    /**
     * Transfer the layout of the given graph to the given associated node.
     *
     * @param node
     *            a compound node
     * @param lgraph
     *            the graph nested in the compound node
     */
    private void graphLayoutToNode(final LNode node, final LGraph lgraph) {
        // Process external ports
        for (LNode childNode : lgraph.getLayerlessNodes()) {
            Object origin = childNode.getProperty(InternalProperties.ORIGIN);
            if (origin instanceof LPort) {
                LPort port = (LPort) origin;
                KVector portPosition =
                        LGraphUtil.getExternalPortPosition(lgraph, childNode, port.getSize().x,
                                port.getSize().y);
                port.getPosition().x = portPosition.x;
                port.getPosition().y = portPosition.y;
                port.setSide(childNode.getProperty(InternalProperties.EXT_PORT_SIDE));
            }
        }

        // Setup the parent node
        KVector actualGraphSize = lgraph.getActualSize();
        if (lgraph.getProperty(InternalProperties.GRAPH_PROPERTIES).contains(
                GraphProperties.EXTERNAL_PORTS)) {
            // Ports have positions assigned
            node.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
            node.getGraph().getProperty(InternalProperties.GRAPH_PROPERTIES)
                    .add(GraphProperties.NON_FREE_PORTS);
            LGraphUtil.resizeNode(node, actualGraphSize, false, true);
        } else {
            // Ports have not been positioned yet - leave this for next layouter
            LGraphUtil.resizeNode(node, actualGraphSize, true, true);
        }
    }

    private boolean isNested(final LGraph graph) {
        return graph.getParentNode() != null;
    }

    // //////////////////////////////////////////////////////////////////////////////
    // Graph Postprocessing (Size and External Ports)

    /**
     * Sets the size of the given graph such that size constraints are adhered to. Furthermore, the padding is
     * added to the graph size and the graph offset. Afterwards, the border spacing property is reset to 0.
     *
     * <p>
     * Major parts of this method are adapted from
     * {@link ElkUtil#resizeNode(org.eclipse.elk.graph.ElkNode, double, double, boolean, boolean)}.
     * </p>
     *
     * <p>
     * Note: This method doesn't care about labels of compound nodes since those labels are not attached to the graph.
     * </p>
     *
     * @param lgraph
     *            the graph to resize.
     */
    private void resizeGraph(final LGraph lgraph) {
        Set<SizeConstraint> sizeConstraint = lgraph.getProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS);
        Set<SizeOptions> sizeOptions = lgraph.getProperty(LayeredOptions.NODE_SIZE_OPTIONS);

        // getActualSize() used to take the border spacing (what is now included in the padding)
        // into account, which is why by this point it had to be cleared since it had already
        // been applied to the offset and the graph size. It currently does not take the padding
        // into account anymore, but if it does, it needs to be cleared again
        KVector calculatedSize = lgraph.getActualSize();
        KVector adjustedSize = new KVector(calculatedSize);

        // calculate the new size
        if (sizeConstraint.contains(SizeConstraint.MINIMUM_SIZE)) {
            KVector minSize = lgraph.getProperty(LayeredOptions.NODE_SIZE_MINIMUM);

            // if minimum width or height are not set, maybe default to default values
            if (sizeOptions.contains(SizeOptions.DEFAULT_MINIMUM_SIZE)) {
                if (minSize.x <= 0) {
                    minSize.x = ElkUtil.DEFAULT_MIN_WIDTH;
                }

                if (minSize.y <= 0) {
                    minSize.y = ElkUtil.DEFAULT_MIN_HEIGHT;
                }
            }

            // apply new size including border spacing
            adjustedSize.x = Math.max(calculatedSize.x, minSize.x);
            adjustedSize.y = Math.max(calculatedSize.y, minSize.y);
        }

        resizeGraphNoReallyIMeanIt(lgraph, calculatedSize, adjustedSize);
    }

    /**
     * Applies a new effective size to a graph that previously had an old size calculated by the layout algorithm. This
     * method takes care of adjusting content alignments as well as external ports that would be misplaced if the new
     * size is larger than the old one.
     *
     * @param lgraph
     *            the graph to apply the size to.
     * @param oldSize
     *            old size as calculated by the layout algorithm.
     * @param newSize
     *            new size that may be larger than the old one.
     */
    private void resizeGraphNoReallyIMeanIt(final LGraph lgraph, final KVector oldSize,
            final KVector newSize) {

        // obey to specified alignment constraints
        Set<ContentAlignment> contentAlignment = lgraph.getProperty(LayeredOptions.CONTENT_ALIGNMENT);

        // horizontal alignment
        if (newSize.x > oldSize.x) {
            if (contentAlignment.contains(ContentAlignment.H_CENTER)) {
                lgraph.getOffset().x += (newSize.x - oldSize.x) / 2f;
            } else if (contentAlignment.contains(ContentAlignment.H_RIGHT)) {
                lgraph.getOffset().x += newSize.x - oldSize.x;
            }
        }

        // vertical alignment
        if (newSize.y > oldSize.y) {
            if (contentAlignment.contains(ContentAlignment.V_CENTER)) {
                lgraph.getOffset().y += (newSize.y - oldSize.y) / 2f;
            } else if (contentAlignment.contains(ContentAlignment.V_BOTTOM)) {
                lgraph.getOffset().y += newSize.y - oldSize.y;
            }
        }

        // correct the position of eastern and southern hierarchical ports, if necessary
        if (lgraph.getProperty(InternalProperties.GRAPH_PROPERTIES).contains(
                GraphProperties.EXTERNAL_PORTS)
                && (newSize.x > oldSize.x || newSize.y > oldSize.y)) {

            // iterate over the graph's nodes, looking for eastern / southern external ports
            // (at this point, the graph's nodes are not divided into layers anymore)
            for (LNode node : lgraph.getLayerlessNodes()) {
                // we're only looking for external port dummies
                if (node.getType() == NodeType.EXTERNAL_PORT) {
                    // check which side the external port is on
                    PortSide extPortSide = node.getProperty(InternalProperties.EXT_PORT_SIDE);
                    if (extPortSide == PortSide.EAST) {
                        node.getPosition().x += newSize.x - oldSize.x;
                    } else if (extPortSide == PortSide.SOUTH) {
                        node.getPosition().y += newSize.y - oldSize.y;
                    }
                }
            }
        }

        // Actually apply the new size
        LPadding padding = lgraph.getPadding();
        lgraph.getSize().x = newSize.x - padding.left - padding.right;
        lgraph.getSize().y = newSize.y - padding.top - padding.bottom;
    }

}
