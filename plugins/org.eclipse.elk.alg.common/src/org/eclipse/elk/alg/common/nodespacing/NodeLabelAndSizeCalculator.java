/*******************************************************************************
 * Copyright (c) 2017, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.nodespacing;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.Cell;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.ContainerArea;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.GridContainerCell;
import org.eclipse.elk.alg.common.nodespacing.internal.NodeContext;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.CellSystemConfigurator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.HorizontalPortPlacementSizeCalculator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.InsidePortLabelCellCreator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.LabelPlacer;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.NodeLabelAndSizeUtilities;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.NodeLabelCellCreator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.NodeSizeCalculator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.PortContextCreator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.PortLabelPlacementCalculator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.PortPlacementCalculator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.VerticalPortPlacementSizeCalculator;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.NodeAdapter;

/**
 * Knows how to calculate the size of a node and how to place its ports. Takes all
 * {@link org.eclipse.elk.core.options.SizeConstraint size constraints} and
 * {@link org.eclipse.elk.core.options.SizeOptions options} into account.
 * 
 * <p>This calculator internally understands the layout of the node as a grid-like system of cells in which different
 * labels and ports are placed. Each cell has a padding and a minimum size, which may or may not be taken into
 * account when calculating the minimum size of the node.</p>
 */
public final class NodeLabelAndSizeCalculator {
    
    /**
     * No instance required.
     */
    private NodeLabelAndSizeCalculator() {
        
    }
    
    
    /**
     * Processes all direct children of the given graph.
     * 
     * @param graph the graph.
     */
    public static void process(final GraphAdapter<?> graph) {
        // Process all of the graph's direct children
        graph.getNodes().forEach(node -> process(graph, node, true, false));
    }
    
    /**
     * Processes the given node which is assumed to be a child of the given graph. Note that this method does not check
     * whether or not this is the case. The worst that can happen, however, is that wrong spacing values are applied
     * during processing.
     * 
     * @param graph
     *            the node's parent graph.
     * @param node
     *            the node to process.
     * @param applyStuff
     *            {@code true} if the node should actually be resized and have its ports and labels positioned,
     *            {@code false} if we should only return the size that would be applied.
     * @param ignoreInsidePortLabels
     *            if {@code true}, we don't place port labels that should be placed inside. This is usually used in
     *            conjunction with {@code applyStuff} by layout algorithms that want to get a lower bound on a
     *            hierarchical node's size, but that handle inside port labels themselves.
     * @return the node's size that was or would be applied.
     */
    public static KVector process(final GraphAdapter<?> graph, final NodeAdapter<?> node, final boolean applyStuff,
            final boolean ignoreInsidePortLabels) {
        
        // Note that, upon Miro's request, each phase of the algorithm was given a code name in the first version of
        // this code. We happily carry on fulfilling this request in this, the second version.

        /* PREPARATORY PREPARATIONS
         * 
         * Create the context objects that hold all of the information relevant to our calculations, including pointers
         * to all the components of the cell system. The different method calls will often just update information in
         * the context object (or nested objects) that subsequent method calls will make use of. Creating the port
         * contexts will also create label cells for each port that has labels.
         */
        NodeContext nodeContext = new NodeContext(graph, node);
        PortContextCreator.createPortContexts(nodeContext, ignoreInsidePortLabels);
        
        
        /* PHASE 1: WONDEROUS WATERFOWL
         *          Setup All Cells
         * 
         * Create all the label cells that will hold node labels, as well as the cell containers that will hold them,
         * for both inside and outside node labels. Also, assign node labels to the relevant label cells. If port
         * labels are to be placed on the inside, setup the inside port label cells with the appropriate paddings, and
         * set the width of the eastern and western ones to the maximum width of the labels they will contain. We can't
         * do that for the northern and southern cells yet because port label placement is more complicated there.
         */
        boolean horizontalLayoutMode = true; 
        // If the graph is null, no layout direction is specified, or the layout direction is set to undefined, 
        // use horizontal layout mode (which yields vertically stacked labels).
        if (graph != null && graph.hasProperty(CoreOptions.DIRECTION)) {
            final Direction layoutDirection = graph.getProperty(CoreOptions.DIRECTION);
            horizontalLayoutMode = layoutDirection == Direction.UNDEFINED || layoutDirection.isHorizontal();  
        }
        NodeLabelCellCreator.createNodeLabelCells(nodeContext, false, horizontalLayoutMode);
        InsidePortLabelCellCreator.createInsidePortLabelCells(nodeContext);
        
        
        /* PHASE 2: DEFECTIVE DUCK
         *          Setup Client Area Space and Node Cell Padding
         * 
         * Apply the minimum client area size to the central grid container cell. Also, reserve space for ports that
         * extend inside the node due to negative port border offsets by setting up appropriate paddings on the main
         * node container cell.
         */
        NodeLabelAndSizeUtilities.setupMinimumClientAreaSize(nodeContext);
        NodeLabelAndSizeUtilities.setupNodePaddingForPortsWithOffset(nodeContext);
        
        
        /* PHASE 3: SALVAGEABLE SWAN
         *          Minimum Space Required to Place Ports
         * 
         * It is now time to find out how large the node needs to be if all ports are to be placed in a way that
         * satisfies all spacing constraints. This may or may not include the labels of ports. We remember these
         * information by setting the minimum width of north / south inside port label cells and the minimum height of
         * east / west inside port label cells. Since the east / west cells are surrounded by the north / south cells,
         * their height may be updated later once we know how high the north / south cells will be.
         */
        HorizontalPortPlacementSizeCalculator.calculateHorizontalPortPlacementSize(nodeContext);
        VerticalPortPlacementSizeCalculator.calculateVerticalPortPlacementSize(nodeContext);
        
        
        /* PHASE 4: DAMNABLE DUCKLING
         *          Setup Cell System Size Contribution Flags
         * 
         * Depending on the size constraints, the different cells may contribute to the height or to the width of the
         * node. In this phase, we setup the size contribution flags according to the size constraints. This lays the
         * groundwork for letting the cell system calculate stuff.
         */
        CellSystemConfigurator.configureCellSystemSizeContributions(nodeContext);
        
        
        /* PHASE 5: DUCK AND COVER
         *          Set Node Width and Place Horizontal Ports
         * 
         * We can now set the node's width and place the ports (and port labels) along the horizontal sides. Since we
         * have no idea how high the node is going to be yet, we place southern ports with the assumption that the
         * node has a height of 0. We will later have to offset those ports by the node's height. Setting the node
         * width has the side effect of computing a horizontal layout for the cell system.
         */
        NodeSizeCalculator.setNodeWidth(nodeContext);
        
        PortPlacementCalculator.placeHorizontalPorts(nodeContext);
        PortLabelPlacementCalculator.placeHorizontalPortLabels(nodeContext);
        
        
        /* PHASE 6: GIGANTIC GOOSE
         *          Set Node Height and Place Vertical Ports
         * 
         * We can now calculate the node's height and place the ports (and port labels) along the vertical sides. Also,
         * since we now know the node's height, we can finally correct the southern port positions. Before we can do
         * all that, however, we might need to update the height and padding of the eastern and western inside port
         * label cells to be sure that free ports are positioned properly.
         * 
         * Note that if we are to not apply stuff, we're done once we know the node's height. Which is why we stop at
         * that point.
         */
        CellSystemConfigurator.updateVerticalInsidePortLabelCellPadding(nodeContext);
        NodeSizeCalculator.setNodeHeight(nodeContext);
        
        if (!applyStuff) {
            return nodeContext.nodeSize;
        }
        
        NodeLabelAndSizeUtilities.offsetSouthernPortsByNodeSize(nodeContext);
        
        PortPlacementCalculator.placeVerticalPorts(nodeContext);
        PortLabelPlacementCalculator.placeVerticalPortLabels(nodeContext);
        
        
        /* PHASE 7: THANKSGIVING
         *          Place Labels and Apply Stuff
         * 
         * Since we now have the node's final size, we can now calculate the positions of the containers for outer node
         * labels and place all inner and outer node labels. Also, the port label cells have positions assigned to them
         * and can be told to apply positions to their labels. Finally, we can apply the node's size and all of the
         * port positions we have calculated.
         */
        LabelPlacer.placeLabels(nodeContext);
        NodeLabelAndSizeUtilities.setNodePadding(nodeContext);
        NodeLabelAndSizeUtilities.applyStuff(nodeContext);
        
        // Return the size
        return nodeContext.nodeSize;
    }
    
    /**
     * Computes the padding required to place inside non-center node labels. This can be used to reserve space around
     * a hierarchical node while laying out its content. The layout direction is interesting here: usually, one would
     * probably either use the graph's layout direction, or the layout direction used to lay out the node's children.
     * However, the way the node's labels are placed may not have anything to do with either direction.
     * 
     * @param graph
     *            the node's parent graph which may specify spacings inherited by the node.
     * @param node
     *            the node to process.
     * @param layoutDirection
     *            the layout direction to assume when computing the paddings.
     * @return the padding required for inside node labels.
     */
    public static ElkPadding computeInsideNodeLabelPadding(final GraphAdapter<?> graph, final NodeAdapter<?> node,
            final Direction layoutDirection) {
        
        // Create a node context and fill it with all the inside node labels
        NodeContext nodeContext = new NodeContext(null, node);
        NodeLabelCellCreator.createNodeLabelCells(nodeContext, true, !layoutDirection.isVertical());
        
        GridContainerCell labelCellContainer = nodeContext.insideNodeLabelContainer;
        ElkPadding padding = new ElkPadding();
        
        // Top
        for (ContainerArea col : ContainerArea.values()) {
            Cell labelCell = labelCellContainer.getCell(ContainerArea.BEGIN, col);
            if (labelCell != null) {
                padding.top = Math.max(padding.top, labelCell.getMinimumHeight());
            }
        }

        // Bottom
        for (ContainerArea col : ContainerArea.values()) {
            Cell labelCell = labelCellContainer.getCell(ContainerArea.END, col);
            if (labelCell != null) {
                padding.bottom = Math.max(padding.bottom, labelCell.getMinimumHeight());
            }
        }

        // Left
        for (ContainerArea row : ContainerArea.values()) {
            Cell labelCell = labelCellContainer.getCell(row, ContainerArea.BEGIN);
            if (labelCell != null) {
                padding.left = Math.max(padding.left, labelCell.getMinimumWidth());
            }
        }

        // Right
        for (ContainerArea row : ContainerArea.values()) {
            Cell labelCell = labelCellContainer.getCell(row, ContainerArea.END);
            if (labelCell != null) {
                padding.right = Math.max(padding.right, labelCell.getMinimumWidth());
            }
        }
        
        // Apply insets and gap where necessary
        if (padding.top > 0) {
            padding.top += labelCellContainer.getPadding().top;
            padding.top += labelCellContainer.getGap();
        }

        if (padding.bottom > 0) {
            padding.bottom += labelCellContainer.getPadding().bottom;
            padding.bottom += labelCellContainer.getGap();
        }

        if (padding.left > 0) {
            padding.left += labelCellContainer.getPadding().left;
            padding.left += labelCellContainer.getGap();
        }

        if (padding.right > 0) {
            padding.right += labelCellContainer.getPadding().right;
            padding.right += labelCellContainer.getGap();
        }
        
        return padding;
    }
    
}
