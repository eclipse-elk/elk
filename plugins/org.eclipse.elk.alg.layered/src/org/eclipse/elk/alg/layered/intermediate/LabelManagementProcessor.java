/*******************************************************************************
 * Copyright (c) 2010, 2017 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.labels.ILabelManager;
import org.eclipse.elk.core.labels.LabelManagementOptions;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Invokes a potentially registered label manager on labels. This processor has two modes of operation. The first mode
 * seeks out node labels, port labels, and edge labels. This mode is used for the first label management run before
 * computing node sizes. The second mode manages center edge labels, which requires node sizes to be fixed.
 * 
 * <dl>
 * <dt>Precondition:</dt>
 * <dd>The graph is layered.</dd>
 * <dd>Label dummy nodes are placed in their final layers.</dd>
 * <dt>Postcondition:</dt>
 * <dd>Possible label text modifications to end labels, port labels, and node labels.</dd>
 * <dt>Slots:</dt>
 * <dd>Before phase 4.</dd>
 * <dt>Same-slot dependencies for mode 1:</dt>
 * <dd>None.</dd>
 * <dt>Same-slot dependencies for mode 2:</dt>
 * <dd>{@link LabelDummySwitcher}</dd>
 * </dl>
 */
public final class LabelManagementProcessor implements ILayoutProcessor<LGraph> {

    /** Minimum width for shortened port labels. */
    private static final double MIN_WIDTH_PORT_LABELS = 20;
    /** Minimum width for shortened node labels. */
    private static final double MIN_WIDTH_NODE_LABELS = 40;
    /** Minimum width for shortened edge labels. */
    public static final double MIN_WIDTH_EDGE_LABELS = 60;

    /** Whether we should only manage edge center labels. */
    private boolean centerLabels;

    /**
     * Creates a new instance.
     * 
     * @param centerLabels
     *            {@code true} if the new instance should operate only on edge center labels, {@code false} if it should
     *            operate on everything except edge center labels.
     */
    public LabelManagementProcessor(final boolean centerLabels) {
        this.centerLabels = centerLabels;
    }

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Label management", 1);

        // If there is no label manager configured for this graph, this processor should not be part of the algorithm.
        // However, better be sure.
        ILabelManager labelManager = layeredGraph.getProperty(LabelManagementOptions.LABEL_MANAGER);
        if (labelManager != null) {
            double edgeLabelSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_LABEL).doubleValue();
            double labelLabelSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_LABEL_LABEL).doubleValue();

            if (centerLabels) {
                manageCenterLabels(layeredGraph, labelManager, edgeLabelSpacing, labelLabelSpacing);
            } else {
                manageNonCenterLabels(layeredGraph, labelManager, labelLabelSpacing);
            }
        }

        monitor.done();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Everything Except Center Edge Labels

    /**
     * Calls label management on all labels that are not edge center labels.
     */
    private void manageNonCenterLabels(final LGraph lGraph, final ILabelManager labelManager,
            final double labelLabelSpacing) {

        boolean verticalLayout = lGraph.getProperty(LayeredOptions.DIRECTION).isVertical();

        // Iterate over the layers
        for (Layer layer : lGraph) {
            // Apply label management to node and port labels
            for (LNode node : layer) {
                if (node.getType() == NodeType.NORMAL) {
                    // Handle node labels
                    doManageLabels(labelManager, node.getLabels(), MIN_WIDTH_NODE_LABELS, labelLabelSpacing,
                            verticalLayout);

                    // Handle ports
                    List<LPort> ports = node.getPorts();
                    for (LPort port : ports) {
                        doManageLabels(labelManager, port.getLabels(), MIN_WIDTH_PORT_LABELS, labelLabelSpacing,
                                verticalLayout);
                    }

                    // Handle attached comments
                    if (node.hasProperty(InternalProperties.TOP_COMMENTS)) {
                        doManageAttachedCommentLabels(labelManager, node.getProperty(InternalProperties.TOP_COMMENTS),
                                MIN_WIDTH_NODE_LABELS, verticalLayout);
                    }

                    if (node.hasProperty(InternalProperties.BOTTOM_COMMENTS)) {
                        doManageAttachedCommentLabels(labelManager,
                                node.getProperty(InternalProperties.BOTTOM_COMMENTS), MIN_WIDTH_NODE_LABELS,
                                verticalLayout);
                    }
                }

                // Edges can have edge Unand tail labels that need to be managed as well (note that at this
                // point, only head and tail labels remain)
                for (LEdge edge : node.getOutgoingEdges()) {
                    doManageLabels(labelManager, edge.getLabels(), MIN_WIDTH_EDGE_LABELS, 0, verticalLayout);
                }
            }
        }
    }

    /**
     * Manages labels of comments attached to a node and thus currently not part of the graph structure.
     */
    private void doManageAttachedCommentLabels(final ILabelManager labelManager, final List<LNode> commentNodes,
            final double minWidthNodeLabels, final boolean verticalLayout) {

        for (LNode commentNode : commentNodes) {
            // We only do stuff if the node actually does have labels
            if (commentNode.getLabels().isEmpty()) {
                continue;
            }

            doManageLabels(labelManager, commentNode.getLabels(), minWidthNodeLabels, 0, verticalLayout);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Center Edge Labels

    /**
     * Calls label management on all edge center labels.
     */
    private void manageCenterLabels(final LGraph lGraph, final ILabelManager labelManager,
            final double edgeLabelSpacing, final double labelLabelSpacing) {

        boolean verticalLayout = lGraph.getProperty(LayeredOptions.DIRECTION).isVertical();

        // Iterate over the layers and find label dummy nodes
        for (Layer layer : lGraph) {
            // The maximum width is used as the target width for center edge labels
            double maxWidth = Math.max(MIN_WIDTH_EDGE_LABELS, LGraphUtil.findMaxNonDummyNodeWidth(layer, false));

            for (LNode node : layer) {
                if (node.getType() == NodeType.LABEL) {
                    LEdge edge = node.getConnectedEdges().iterator().next();
                    double edgeThickness = edge.getProperty(LayeredOptions.EDGE_THICKNESS).doubleValue();

                    // The list of labels should never be empty (otherwise the label dummy wouldn't have been created
                    // in the first place)
                    Iterable<LLabel> labels = node.getProperty(InternalProperties.REPRESENTED_LABELS);
                    KVector spaceRequiredForLabels =
                            doManageLabels(labelManager, labels, maxWidth, labelLabelSpacing, verticalLayout);

                    // Apply the space required for labels to the dummy node (we don't bother with the ports here since
                    // they will be meddled with later by the LabelSideSelector anyway)
                    node.getSize().x = spaceRequiredForLabels.x;
                    node.getSize().y = spaceRequiredForLabels.y + edgeThickness + edgeLabelSpacing;
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Actual Label Management

    /**
     * Manage all the labels according to the labelManager and change the labels size.
     * 
     * @param labelManager
     *            the label manager used to manage labels.
     * @param labels
     *            the labels to be passed to the label manager.
     * @param targetWidth
     *            the target width.
     * @param labelLabelSpacing
     *            spacing between two adjacent labels
     * @param verticalLayout
     *            {@code true} if the graph is laid out downwards or upwards. In these cases, the new label dimensions
     *            returned by the label manager need to be turned 90 degrees.
     * @return size required to place the labels. Does not include the edge-label space yet.
     */
    public static KVector doManageLabels(final ILabelManager labelManager, final Iterable<LLabel> labels,
            final double targetWidth, final double labelLabelSpacing, final boolean verticalLayout) {

        KVector requiredLabelSpace = new KVector();

        // Manage all of them labels, if there are any
        if (labels.iterator().hasNext()) {
            for (LLabel label : labels) {
                KVector labelSize = label.getSize();

                // If the label has an origin, call the label size modifier
                Object origin = label.getProperty(InternalProperties.ORIGIN);
                if (origin != null) {
                    KVector newSize = labelManager.manageLabelSize(origin, targetWidth);

                    if (newSize != null) {
                        if (verticalLayout) {
                            // Our labels are turned 90 degrees
                            labelSize.x = newSize.y;
                            labelSize.y = newSize.x;
                        } else {
                            labelSize.x = newSize.x;
                            labelSize.y = newSize.y;
                        }
                    }
                }

                // Update required label space
                if (verticalLayout) {
                    requiredLabelSpace.x += labelLabelSpacing + label.getSize().x;
                    requiredLabelSpace.y = Math.max(requiredLabelSpace.y, labelSize.y);
                } else {
                    requiredLabelSpace.x = Math.max(requiredLabelSpace.x, labelSize.x);
                    requiredLabelSpace.y += labelLabelSpacing + label.getSize().y;
                }
            }

            // There is one label-label spacing too much
            if (verticalLayout) {
                requiredLabelSpace.x -= labelLabelSpacing;
            } else {
                requiredLabelSpace.y -= labelLabelSpacing;
            }
        }

        return requiredLabelSpace;
    }

}
