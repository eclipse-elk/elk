/*******************************************************************************
 * Copyright (c) 2012, 2020 Kiel University and others.
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

import java.util.EnumSet;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.HorizontalLabelAlignment;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.LabelCell;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.VerticalLabelAlignment;
import org.eclipse.elk.alg.common.overlaps.RectangleStripOverlapRemover;
import org.eclipse.elk.alg.common.overlaps.RectangleStripOverlapRemover.OverlapRemovalDirection;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphAdapters;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LMargin;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.core.options.LabelSide;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * <p>Puts all end labels into {@link LabelCell label cells} stored for each node using the
 * {@link InternalProperties#END_LABELS} property. The label cells are properly positioned relative to their node's
 * top left corner, and thus will need post-processing once node coordinates are known. Also updates node margins to
 * include end labels.</p>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a layered graph</dd>
 *     <dd>labels are marked with a placement side</dd>
 *     <dd>node sizes and port positions are fixed</dd>
 *     <dd>port lists are sorted according to {@link PortListSorter#PortComparator}</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>each node with end labels has its {@link InternalProperties#END_LABELS} property set accordingly</dd>
 *     <dd>label cells are large enough to hold their labels and are positioned correctly relative to their node</dd>
 *     <dd>node margins include space required to place edge end labels.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 4.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link LabelAndNodeSizeProcessor}</dd>
 *     <dd>{@link InnermostNodeMarginCalculator}</dd>
 *     <dd>{@link LabelSideSelector}</dd>
 * </dl>
 * 
 * @see EndLabelPostprocessor
 */
public final class EndLabelPreprocessor implements ILayoutProcessor<LGraph> {
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("End label pre-processing", 1);
        
        double edgeLabelSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_LABEL);
        double labelLabelSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_LABEL_LABEL);
        boolean verticalLayout = layeredGraph.getProperty(LayeredOptions.DIRECTION).isVertical();
        
        // We iterate over each node and place the end labels of its incident edges
        layeredGraph.getLayers().stream()
                .flatMap(layer -> layer.getNodes().stream())
                .forEach(node -> processNode(node, edgeLabelSpacing, labelLabelSpacing, verticalLayout));
        
        monitor.done();
    }
    
    private void processNode(final LNode node, final double edgeLabelSpacing, final double labelLabelSpacing,
            final boolean verticalLayout) {

        // Iterate over all ports and collect their labels in label cells
        int portCount = node.getPorts().size(); 
        LabelCell[] portLabelCells = new LabelCell[portCount];
        
        for (int portIndex = 0; portIndex < portCount; portIndex++) {
            LPort port = node.getPorts().get(portIndex);
            port.id = portIndex;
            
            portLabelCells[portIndex] = createConfiguredLabelCell(
                    gatherLabels(port), labelLabelSpacing, verticalLayout);
        }
        
        // Actually go off and place them labels!
        placeLabels(node, portLabelCells, labelLabelSpacing, edgeLabelSpacing, verticalLayout);
        
        // Turn the array into a map and save that in the node
        Map<LPort, LabelCell> portToLabelCellMap = new HashMap<>();
        for (int index = 0; index < portLabelCells.length; index++) {
            if (portLabelCells[index] != null) {
                portToLabelCellMap.put(node.getPorts().get(index), portLabelCells[index]);
            }
        }
        
        if (!portToLabelCellMap.isEmpty()) {
            node.setProperty(InternalProperties.END_LABELS, portToLabelCellMap);
            
            // Update the node's margins
            updateNodeMargins(node, portLabelCells);
        }
    }
    
    /**
     * Creates label cell for the given port with the given labels, if any. If there are no labels, this method returns
     * {@code null}. The label cell will still have to have its alignment set up, but its size is already set to the
     * size required to place its labels.
     */
    private LabelCell createConfiguredLabelCell(final List<LLabel> labels, final double labelLabelSpacing,
            final boolean verticalLayout) {
        
        if (labels == null || labels.isEmpty()) {
            return null;
        }
        
        // Create the new label cell and setup its alignments depending on the port's side
        LabelCell labelCell = new LabelCell(labelLabelSpacing, !verticalLayout);
        
        for (LLabel label : labels) {
            labelCell.addLabel(LGraphAdapters.adapt(label));
        }
        
        // Setup the label cell's size
        ElkRectangle labelCellRect = labelCell.getCellRectangle();
        labelCellRect.height = labelCell.getMinimumHeight();
        labelCellRect.width = labelCell.getMinimumWidth();
        
        return labelCell;
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Label Gathering
    
    /** Special value to indicate that there are no edges incident to a port. */
    private static final double NO_INCIDENT_EDGE_THICKNESS = -1;

    /**
     * Returns a list that contains all end labels to be placed at the given port. If there are no such labels, the
     * list will be empty. If there are not even any edges connected to the port, the return value will be {@code null}.
     * 
     *  <p>his method also takes care of north / south ports whose original edges have been rerouted to a north / south
     *  port dummy. The maximum thickness of any edge connected to the port is saved
     *  {@link InternalProperties#MAX_EDGE_THICKNESS as a property} if it has labels.</p>
     * 
     * <p>This method is also used by {@link LabelSideSelector}.</p>
     */
    static List<LLabel> gatherLabels(final LPort port) {
        List<LLabel> labels = Lists.newArrayList();
        
        // Gather labels of the port itself
        double maxEdgeThickness = gatherLabels(port, labels);
        
        // If it has a dummy associated with it, we need to go through the dummy's ports and process those that were
        // created for the current port (see the NorthSouthPortPreprocessor)
        LNode dummyNode = port.getProperty(InternalProperties.PORT_DUMMY);
        if (dummyNode != null) {
            for (LPort dummyPort : dummyNode.getPorts()) {
                if (dummyPort.getProperty(InternalProperties.ORIGIN) == port) {
                    maxEdgeThickness = Math.max(
                            maxEdgeThickness,
                            gatherLabels(dummyPort, labels));
                }
            }
        }
        
        // Only save the maximum edge thickness if we'll be interested in it later
        if (!labels.isEmpty()) {
            port.setProperty(InternalProperties.MAX_EDGE_THICKNESS, maxEdgeThickness);
        }
        
        return maxEdgeThickness != NO_INCIDENT_EDGE_THICKNESS ? labels : null;
    }
    
    /**
     * Puts all relevant end labels of edges connected to the given port into the given list. Returns the maximum edge
     * thickness of any incident edge or {@link #NO_INCIDENT_EDGE_THICKNESS} if there is no incident edge.
     */
    private static double gatherLabels(final LPort port, final List<LLabel> targetList) {
        double maxEdgeThickness = -1;
        
        List<LLabel> labels = new LinkedList<>();
        
        for (LEdge incidentEdge : port.getConnectedEdges()) {
            maxEdgeThickness = Math.max(maxEdgeThickness, incidentEdge.getProperty(LayeredOptions.EDGE_THICKNESS));
            
            if (incidentEdge.getSource() == port) {
                // It's an outgoing edge; all tail labels belong to this port
                incidentEdge.getLabels().stream()
                        .filter(label -> label.getProperty(
                                LayeredOptions.EDGE_LABELS_PLACEMENT) == EdgeLabelPlacement.TAIL)
                        .forEach(label -> labels.add(label));
            } else {
                // It's an incoming edge; all head labels belong to this port
                incidentEdge.getLabels().stream()
                        .filter(label -> label.getProperty(
                                LayeredOptions.EDGE_LABELS_PLACEMENT) == EdgeLabelPlacement.HEAD)
                        .forEach(label -> labels.add(label));
            }
            
            // Remember the edge each label came from
            for (LLabel label : labels) {
                if (!label.hasProperty(InternalProperties.END_LABEL_EDGE)) {
                    label.setProperty(InternalProperties.END_LABEL_EDGE, incidentEdge);
                }
            }
            
            targetList.addAll(labels);
            labels.clear();
        }
        
        return maxEdgeThickness;
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Label Placement
    
    /**
     * Places end labels of all of the node's ports.
     */
    private void placeLabels(final LNode node, final LabelCell[] portLabelCells, final double labelLabelSpacing,
            final double edgeLabelSpacing, final boolean verticalLayout) {
        
        // First, place them as we usually would. This step can result in overlaps for northern / southern labels
        // (for horizontal layout directions) or for eastern / western labels (for vertical layout directions).
        EnumSet<PortSide> portSidesWithLabels = EnumSet.noneOf(PortSide.class);
        for (LPort port : node.getPorts()) {
            if (portLabelCells[port.id] != null) {
                placeLabels(port, portLabelCells[port.id], edgeLabelSpacing);
                portSidesWithLabels.add(port.getSide());
            }
        }
        
        // If there are ports on the problematic sides, go ahead and remove overlaps between them
        if (verticalLayout) {
            removeLabelOverlaps(node, portLabelCells, PortSide.EAST, 2 * labelLabelSpacing, edgeLabelSpacing);
            removeLabelOverlaps(node, portLabelCells, PortSide.WEST, 2 * labelLabelSpacing, edgeLabelSpacing);
        } else {
            removeLabelOverlaps(node, portLabelCells, PortSide.NORTH, 2 * labelLabelSpacing, edgeLabelSpacing);
            removeLabelOverlaps(node, portLabelCells, PortSide.SOUTH, 2 * labelLabelSpacing, edgeLabelSpacing);
        }
    }

    /**
     * Places the edge end labels that are to be placed near the given port.
     */
    private void placeLabels(final LPort port, final LabelCell labelCell, final double edgeLabelSpacing) {
        // Some necessary position information
        ElkRectangle labelCellRect = labelCell.getCellRectangle();
        KVector nodeSize = port.getNode().getSize();
        LMargin nodeMargin = port.getNode().getMargin();
        KVector portPos = port.getPosition();
        KVector portAnchor = KVector.sum(portPos, port.getAnchor());
        
        // Calculate cell position depending on port side
        switch (port.getSide()) {
        case NORTH:
            labelCell.setVerticalAlignment(VerticalLabelAlignment.BOTTOM);
            labelCellRect.y = -nodeMargin.top
                    - edgeLabelSpacing
                    - labelCellRect.height;
            
            if (getLabelSide(labelCell) == LabelSide.ABOVE) {
                labelCell.setHorizontalAlignment(HorizontalLabelAlignment.RIGHT);
                labelCellRect.x = portAnchor.x
                        - maxEdgeThickness(port)
                        - edgeLabelSpacing
                        - labelCellRect.width;
            } else {
                labelCell.setHorizontalAlignment(HorizontalLabelAlignment.LEFT);
                labelCellRect.x = portAnchor.x
                        + maxEdgeThickness(port)
                        + edgeLabelSpacing;
            }
            break;
            
        case EAST:
            labelCell.setHorizontalAlignment(HorizontalLabelAlignment.LEFT);
            labelCellRect.x = nodeSize.x
                    + nodeMargin.right
                    + edgeLabelSpacing;
            
            if (getLabelSide(labelCell) == LabelSide.ABOVE) {
                labelCell.setVerticalAlignment(VerticalLabelAlignment.BOTTOM);
                labelCellRect.y = portAnchor.y
                        - maxEdgeThickness(port)
                        - edgeLabelSpacing
                        - labelCellRect.height;
            } else {
                labelCell.setVerticalAlignment(VerticalLabelAlignment.TOP);
                labelCellRect.y = portAnchor.y
                        + maxEdgeThickness(port)
                        + edgeLabelSpacing;
            }
            break;
            
        case SOUTH:
            labelCell.setVerticalAlignment(VerticalLabelAlignment.TOP);
            labelCellRect.y = nodeSize.y
                    + nodeMargin.bottom
                    + edgeLabelSpacing;
            
            if (getLabelSide(labelCell) == LabelSide.ABOVE) {
                labelCell.setHorizontalAlignment(HorizontalLabelAlignment.RIGHT);
                labelCellRect.x = portAnchor.x
                        - maxEdgeThickness(port)
                        - edgeLabelSpacing
                        - labelCellRect.width;
            } else {
                labelCell.setHorizontalAlignment(HorizontalLabelAlignment.LEFT);
                labelCellRect.x = portAnchor.x
                        + maxEdgeThickness(port)
                        + edgeLabelSpacing;
            }
            break;
            
        case WEST:
            labelCell.setHorizontalAlignment(HorizontalLabelAlignment.RIGHT);
            labelCellRect.x = -nodeMargin.left
                    - edgeLabelSpacing
                    - labelCellRect.width;
            
            if (getLabelSide(labelCell) == LabelSide.ABOVE) {
                labelCell.setVerticalAlignment(VerticalLabelAlignment.BOTTOM);
                labelCellRect.y = portAnchor.y
                        - maxEdgeThickness(port)
                        - edgeLabelSpacing
                        - labelCellRect.height;
            } else {
                labelCell.setVerticalAlignment(VerticalLabelAlignment.TOP);
                labelCellRect.y = portAnchor.y
                        + maxEdgeThickness(port)
                        + edgeLabelSpacing;
            }
            break;
        }
    }

    /**
     * Calls the rectangle overlap removal code to remove overlaps between end labels of edges connected to ports on
     * the given sides.
     */
    private void removeLabelOverlaps(final LNode node, final LabelCell[] portLabelCells, final PortSide portSide,
            final double labelLabelSpacing, final double edgeLabelSpacing) {
        
        RectangleStripOverlapRemover overlapRemover = RectangleStripOverlapRemover
                .createForDirection(portSideToOverlapRemovalDirection(portSide))
                .withGap(labelLabelSpacing, labelLabelSpacing)
                .withStartCoordinate(calculateOverlapStartCoordinate(node, portSide, edgeLabelSpacing));
        
        // Gather the rectangles
        for (LPort port : node.getPorts(portSide)) {
            if (portLabelCells[port.id] != null) {
                ElkRectangle labelCellRect = portLabelCells[port.id].getCellRectangle();
                assert labelCellRect.height > 0 && labelCellRect.width > 0;
                
                overlapRemover.addRectangle(labelCellRect);
            }
        }
        
        // Remove overlaps. Since we have stuffed the original label cell rectangles into this thing, we won't evern
        // have to apply or post-process anything. The marvels of modern technology.
        overlapRemover.removeOverlaps();
    }
    
    /**
     * Calculates the start coordinate to use for overlap removal on the given port side.
     */
    private double calculateOverlapStartCoordinate(final LNode node, final PortSide portSide,
            final double edgeLabelSpacing) {
        
        KVector nodeSize = node.getSize();
        LMargin nodeMargin = node.getMargin();
        
        switch (portSide) {
        case NORTH:
            return -nodeMargin.top - edgeLabelSpacing;
        case SOUTH:
            return nodeSize.y + nodeMargin.bottom + edgeLabelSpacing;
        case EAST:
            return nodeSize.x + nodeMargin.right + edgeLabelSpacing;
        case WEST:
            return -nodeMargin.left - edgeLabelSpacing;
        default:
            // Should never happen
            assert false;
            return 0;
        }
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Node Margins
    
    /**
     * Updates the node's margins to account for its end labels.
     */
    private void updateNodeMargins(final LNode node, final LabelCell[] labelCells) {
        LMargin nodeMargin = node.getMargin();
        KVector nodeSize = node.getSize();
        
        // Calculate the rectangle that describes the node's current margin
        ElkRectangle nodeMarginRectangle = new ElkRectangle(
                -nodeMargin.left,
                -nodeMargin.top,
                nodeMargin.left + nodeSize.x + nodeMargin.right,
                nodeMargin.top + nodeSize.y + nodeMargin.bottom);
        
        // Union the rectangle with each rectangle that describes a label cell
        for (LabelCell labelCell : labelCells) {
            if (labelCell != null) {
                nodeMarginRectangle.union(labelCell.getCellRectangle());
            }
        }
        
        // Reapply the new rectangle to the margin
        nodeMargin.left = -nodeMarginRectangle.x;
        nodeMargin.top = -nodeMarginRectangle.y;
        nodeMargin.right = nodeMarginRectangle.width - nodeMargin.left - nodeSize.x;
        nodeMargin.bottom = nodeMarginRectangle.height - nodeMargin.top - nodeSize.y;
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utility Methods
    
    /**
     * Retrieve the side of the edge the labels of the given cell should be placed at. Different labels may actually
     * have different label sides configured for them, but we simply put them all into one place. That that!
     */
    private LabelSide getLabelSide(final LabelCell labelCell) {
        assert labelCell != null;
        assert labelCell.hasLabels();
        
        return labelCell.getLabels().get(0).getProperty(InternalProperties.LABEL_SIDE);
    }
    
    /**
     * Returns the maximum thickness of all edges incident to the port.
     */
    private double maxEdgeThickness(final LPort port) {
        return port.getProperty(InternalProperties.MAX_EDGE_THICKNESS);
    }
    
    /**
     * Returns the overlap removal direction appropriate for the given port side.
     */
    private OverlapRemovalDirection portSideToOverlapRemovalDirection(final PortSide portSide) {
        switch (portSide) {
        case NORTH:
            return OverlapRemovalDirection.UP;
        case SOUTH:
            return OverlapRemovalDirection.DOWN;
        case EAST:
            return OverlapRemovalDirection.RIGHT;
        case WEST:
            return OverlapRemovalDirection.LEFT;
        default:
            // Shouldn't happen
            assert false;
            return null;
        }
    }

}
