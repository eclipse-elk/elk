/*******************************************************************************
 * Copyright (c) 2012, 2017 Kiel University and others.
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

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import org.eclipse.elk.alg.common.nodespacing.NodeLabelAndSizeCalculator;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.CenterEdgeLabelPlacementStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.Alignment;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;

import com.google.common.collect.Lists;


/**
 * Processor that tries to move label dummy nodes into an "optimal" layer their long edges crosses. This
 * is done by switching the order of long edge dummies and label dummies. There are
 * {@link CenterEdgeLabelPlacementStrategy different strategies} available to determine the target layer.
 * 
 * <p>
 * If this was the only thing we did we could end up in situations where multiple edges forming a
 * hyperedge are merged such that it's not clear anymore which edge label belongs to which edge:
 * </p>
 * <pre>
 *       An edge label
 *    -------------------+-----------+-------------------------- - - -
 *                       |           |
 *    -------------------+           +-------------------------- - - -
 *                                       Another edge label
 * </pre>
 * <p>
 * We solve this by making sure that the long edge dummies preceding a label dummy node have their
 * {@link InternalProperties#LONG_EDGE_TARGET} property set to {@code null} (we more or less view the
 * label dummy as the new target of those long edge dummies). The same is true for the
 * {@link InternalProperties#LONG_EDGE_SOURCE} property of long edge dummies succeeding a label dummy.
 * </p>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a properly layered graph</dd>
 *     <dd>center edge labels are represented by center label dummy nodes</dd>
 *     <dd>each label dummy and long edge dummy has exactly one incoming and one outgoing edge</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>center label dummy nodes are the centermost dummies of a long edge</dd>
 *     <dd>the {@link InternalProperties#LONG_EDGE_TARGET} property of long edge dummies preceding center
 *         edge label dummies are set to {@code null} to prevent incorrect hyperedge dummy merging; the
 *         same is true for the {@link InternalProperties#LONG_EDGE_SOURCE} property of long edge dummies
 *         succeeding center edge label dummies.
 *   <dt>Slots:</dt>
 *     <dd>Before phase 4.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link NodeLabelAndSizeCalculator}</dd>
 *     <dd>{@link InnermostNodeMarginCalculator}</dd>
 *     <dd>{@link EndLabelPreprocessor}</dd>
 * </dl>
 * 
 * @see CenterEdgeLabelPlacementStrategy
 */
public final class LabelDummySwitcher implements ILayoutProcessor<LGraph> {
    
    /** A property to communicate with the analyses that can be run on a graph. */
    public static final IProperty<Boolean> INCLUDE_LABEL =
            new Property<>("edgelabelcenterednessanalysis.includelabel", Boolean.FALSE);
    
    /** Width of each layer. */
    private double[] layerWidths = null;
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Label dummy switching", 1);
        
        // The default placement strategy to be assumed if no more specific strategy is
        CenterEdgeLabelPlacementStrategy defaultPlacementStrategy =
                layeredGraph.getProperty(LayeredOptions.EDGE_LABELS_CENTER_LABEL_PLACEMENT_STRATEGY);
        
        // Assign layer IDs
        assignIdsToLayers(layeredGraph);
        
        // Collect label dummy infos by placement strategy
        Map<CenterEdgeLabelPlacementStrategy, List<LabelDummyInfo>> labelDummyInfos =
                gatherLabelDummyInfos(layeredGraph, defaultPlacementStrategy);
        
        // If at least one width-based strategy is active, calculate layers widths. After this point, layer widths are
        // either all zero or equal to each layer's widest non-dummy node
        layerWidths = new double[layeredGraph.getLayers().size()];
        
        for (CenterEdgeLabelPlacementStrategy strategy : CenterEdgeLabelPlacementStrategy.values()) {
            if (strategy.usesLabelSizeInformation() && !labelDummyInfos.get(strategy).isEmpty()) {
                // We have found an active width-based strategy
                calculateLayerWidths(layeredGraph);
                break;
            }
        }
        
        // Work through the non-width-based strategies first. They might change the width of layers, which might
        // influence size-based strategies later on
        for (CenterEdgeLabelPlacementStrategy strategy : CenterEdgeLabelPlacementStrategy.values()) {
            if (!strategy.usesLabelSizeInformation()) {
                processStrategy(labelDummyInfos.get(strategy));
            }
        }
        
        // Now execute size-based strategies
        for (CenterEdgeLabelPlacementStrategy strategy : CenterEdgeLabelPlacementStrategy.values()) {
            if (strategy.usesLabelSizeInformation()) {
                processStrategy(labelDummyInfos.get(strategy));
            }
        }
        
        // Reset data
        layerWidths = null;
        
        monitor.done();
    }
    
    /**
     * Goes through the graph's layers and assigns zero-based IDs.
     */
    private void assignIdsToLayers(final LGraph layeredGraph) {
        int layerIndex = 0;
        for (Layer layer : layeredGraph) {
            layer.id = layerIndex;
            layerIndex++;
        }
    }
    
    /**
     * Gathers the graph's label dummies, turns them into {@link LabelDummyInfo} objects, and returns them by placement
     * strategy.
     */
    private Map<CenterEdgeLabelPlacementStrategy, List<LabelDummyInfo>> gatherLabelDummyInfos(
            final LGraph layeredGraph, final CenterEdgeLabelPlacementStrategy defaultPlacementStrategy) {
        
        // Initialize the map with empty lists
        Map<CenterEdgeLabelPlacementStrategy, List<LabelDummyInfo>> infos =
                new EnumMap<>(CenterEdgeLabelPlacementStrategy.class);
        
        for (CenterEdgeLabelPlacementStrategy strategy : CenterEdgeLabelPlacementStrategy.values()) {
            infos.put(strategy, new ArrayList<LabelDummyInfo>());
        }
        
        // From the graph's nodes, filter out the label dummies and add them to the appropriate lists
        layeredGraph.getLayers().stream()
                .flatMap((layer) -> layer.getNodes().stream())
                .filter(node -> node.getType() == NodeType.LABEL)
                .map(labelDummy -> new LabelDummyInfo(labelDummy, defaultPlacementStrategy))
                .forEach(dummyInfo -> infos.get(dummyInfo.placementStrategy).add(dummyInfo));
        
        return infos;
    }

    /**
     * Fills the layer widths array, indexed by layer ID which was assigned to the layers previously.
     */
    private void calculateLayerWidths(final LGraph layeredGraph) {
        assert layerWidths.length == layeredGraph.getLayers().size();
        
        for (Layer layer : layeredGraph.getLayers()) {
            layerWidths[layer.id] = LGraphUtil.findMaxNonDummyNodeWidth(layer, false);
        }
    }
    
    /**
     * Executes label dummy switching on the given list of label dummy infos. For most strategies, this method simply
     * invokes the strategy on each label dummy. However, if a strategy places all dummies at once, this method
     * completely delegates to that strategy's implementation.
     */
    private void processStrategy(final List<LabelDummyInfo> labelDummyInfos) {
        // There might not be anything to do for this strategy
        if (labelDummyInfos.isEmpty()) {
            return;
        }
        
        // Check if the strategy has a special processing method
        if (labelDummyInfos.get(0).placementStrategy == CenterEdgeLabelPlacementStrategy.SPACE_EFFICIENT_LAYER) {
            computeSpaceEfficientAssignment(labelDummyInfos);
            
        } else {
            // Execute the strategy for each label dummy
            for (LabelDummyInfo labelDummyInfo : labelDummyInfos) {
                switch (labelDummyInfo.placementStrategy) {
                case CENTER_LAYER:
                    assignLayer(labelDummyInfo, findCenterLayerTargetId(labelDummyInfo));
                    break;
                    
                case MEDIAN_LAYER:
                    assignLayer(labelDummyInfo, findMedianLayerTargetId(labelDummyInfo));
                    break;
                    
                case WIDEST_LAYER:
                    assignLayer(labelDummyInfo, findWidestLayerTargetId(labelDummyInfo));
                    break;
                    
                case HEAD_LAYER:
                    setEndLayerNodeAlignment(labelDummyInfo);
                    assignLayer(labelDummyInfo, findEndLayerTargetId(labelDummyInfo, true));
                    break;
                    
                case TAIL_LAYER:
                    setEndLayerNodeAlignment(labelDummyInfo);
                    assignLayer(labelDummyInfo, findEndLayerTargetId(labelDummyInfo, false));
                    break;
                }
                
                updateLongEdgeSourceLabelDummyInfo(labelDummyInfo);
            }
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Widest Layer

    /**
     * Find a layer to place the label dummy into. This method bases its decision on which long edge dummy is in the
     * widest layer.
     */
    private int findWidestLayerTargetId(final LabelDummyInfo labelDummyInfo) {
        // Find the widest layer among those the long edge dummies are placed in
        int widestLayerIndex = labelDummyInfo.leftmostLayerId;
        
        for (int index = widestLayerIndex + 1; index <= labelDummyInfo.rightmostLayerId; index++) {
            if (layerWidths[index] > layerWidths[widestLayerIndex]) {
                widestLayerIndex = index;
            }
        }
        
        // Return the dummy node in the widest layer
        return widestLayerIndex;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Center Layer

    /**
     * Find a layer to place the label dummy into. This method bases its decision on which layer it thinks will be
     * closest to the edge's physical center later.
     */
    private int findCenterLayerTargetId(final LabelDummyInfo labelDummyInfo) {
        // Sum up the widths of all the layers this thing spans
        double[] layerWidthSums = computeLayerWidthSums(labelDummyInfo);
        
        // Find the first layer that exceeds half the width
        double threshold = layerWidthSums[layerWidthSums.length - 1] / 2;
        
        for (int i = 0; i < layerWidthSums.length; i++) {
            if (layerWidthSums[i] >= threshold) {
                return labelDummyInfo.leftmostLayerId + i;
            }
        }
        
        // This should actually not happen
        assert false;
        return labelDummyInfo.leftmostLayerId + labelDummyInfo.leftLongEdgeDummies.size();
    }
    
    /**
     * Computes an array in which entry i refers to the combined width of layers [0, i]. The combined width is
     * estimated based on the currently estimated layer width and a wild guess as to the layer spacing.
     */
    private double[] computeLayerWidthSums(final LabelDummyInfo labelDummyInfo) {
        // The minimum space that we think will be left between 
        LGraph lgraph = labelDummyInfo.labelDummy.getGraph();
        double edgeNodeSpacing = lgraph.getProperty(LayeredOptions.SPACING_EDGE_NODE_BETWEEN_LAYERS) * 2;
        double nodeNodeSpacing = lgraph.getProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS);
        double minSpaceBetweenLayers = Math.max(edgeNodeSpacing, nodeNodeSpacing);
        
        // The array that will hold the accumulated widths
        double[] layerWidthSums = new double[labelDummyInfo.totalDummyCount()];
        
        double currentWidthSum = -minSpaceBetweenLayers;
        int currentIndex = 0;
        
        for (LNode leftDummy : labelDummyInfo.leftLongEdgeDummies) {
            currentWidthSum += layerWidths[leftDummy.getLayer().id] + minSpaceBetweenLayers;
            layerWidthSums[currentIndex++] = currentWidthSum;
        }
        
        currentWidthSum += layerWidths[labelDummyInfo.labelDummy.getLayer().id] + minSpaceBetweenLayers;
        layerWidthSums[currentIndex++] = currentWidthSum;
        
        for (LNode rightDummy : labelDummyInfo.rightLongEdgeDummies) {
            currentWidthSum += layerWidths[rightDummy.getLayer().id] + minSpaceBetweenLayers;
            layerWidthSums[currentIndex++] = currentWidthSum;
        }
        
        return layerWidthSums;
    }

    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Median Layer

    /**
     * Find a layer to place the label dummy into. This method bases its decision on which long edge dummy is in the
     * median layer of the layers the long edge spans.
     */
    private int findMedianLayerTargetId(final LabelDummyInfo labelDummyInfo) {
        // Find the median of the layers spanned by the long edge this label dummy is part of
        int layers = labelDummyInfo.totalDummyCount();
        int lowerMedian = (layers - 1) / 2;
        
        return labelDummyInfo.leftmostLayerId + lowerMedian;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // End Layer

    /**
     * Find a layer to place the label dummy into. This method simply returns the dummy node closest to the original
     * edge's head or tail. It does take into account whether the edge is reversed or not.
     * 
     * @param headLayer
     *            {@code true} if the label is to be placed at the edge's head, {@code false} if it is to be placed at
     *            the edge's tail.
     */
    private int findEndLayerTargetId(final LabelDummyInfo labelDummyInfo, final boolean headLayer) {
        boolean reversed = isPartOfReversedEdge(labelDummyInfo);
        
        if ((headLayer && !reversed) || (!headLayer && reversed)) {
            return labelDummyInfo.rightmostLayerId;
        } else {
            return labelDummyInfo.leftmostLayerId;
        }
    }
    
    /**
     * Sets the alignment property for the given label dummy. Its placement strategy is assumed to be one of the end
     * strategies.
     */
    private void setEndLayerNodeAlignment(final LabelDummyInfo labelDummyInfo) {
        assert labelDummyInfo.placementStrategy == CenterEdgeLabelPlacementStrategy.HEAD_LAYER
                || labelDummyInfo.placementStrategy == CenterEdgeLabelPlacementStrategy.TAIL_LAYER;
        
        boolean isHeadLabel = labelDummyInfo.placementStrategy == CenterEdgeLabelPlacementStrategy.HEAD_LAYER;
        boolean isPartOfReversedEdge = isPartOfReversedEdge(labelDummyInfo);
        
        if ((isHeadLabel && !isPartOfReversedEdge) || (!isHeadLabel && isPartOfReversedEdge)) {
            labelDummyInfo.labelDummy.setProperty(LayeredOptions.ALIGNMENT, Alignment.RIGHT);
        } else {
            labelDummyInfo.labelDummy.setProperty(LayeredOptions.ALIGNMENT, Alignment.LEFT);
        }
    }
    
    /**
     * Checks if the given label dummy node is part of a reversed edge.
     */
    private boolean isPartOfReversedEdge(final LabelDummyInfo labelDummyInfo) {
        assert labelDummyInfo.labelDummy.getType() == NodeType.LABEL;
        assert labelDummyInfo.labelDummy.getIncomingEdges().iterator().hasNext();
        assert labelDummyInfo.labelDummy.getOutgoingEdges().iterator().hasNext();
        
        // Find incoming and outgoing edge
        LEdge incoming = labelDummyInfo.labelDummy.getIncomingEdges().iterator().next();
        LEdge outgoing = labelDummyInfo.labelDummy.getOutgoingEdges().iterator().next();
        
        return incoming.getProperty(InternalProperties.REVERSED) || outgoing.getProperty(InternalProperties.REVERSED);
    }

    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Space Efficient
    
    /* The space-efficient assignment strategy is a heuristic. We first perform all trivial assignments: labels
     * that only have one valid layer available to them, or that can be assigned to a layer larger than the label.
     * We are then left with labels that are too big for their layers. Given a current label, we define a layer's
     * potential width to be the maximum of its current width and the width of all unassigned labels that can be
     * assigned to the layer, except for the current label. We iterate over our labels from biggest to smallest and
     * assign each to the layer with the currently largest potential width.
     */
    
    /**
     * Runs a heuristic that tries to compute the least wide label dummy to layer assignment.
     */
    private void computeSpaceEfficientAssignment(final List<LabelDummyInfo> labelDummyInfos) {
        // We start by assigning all label dummies that only have a single layer to choose from or that can be
        // assigned to a layer large enough for them and removing them from our list
        List<LabelDummyInfo> nonTrivialLabels = performTrivialAssignments(labelDummyInfos);
        if (nonTrivialLabels.isEmpty()) {
            return;
        }
        
        // The remaining labels are not as easy to assign. Sort descendingly by size.
        nonTrivialLabels.sort(
                (info1, info2) -> Double.compare(info2.labelDummy.getSize().x, info1.labelDummy.getSize().x));
        
        int labelCount = nonTrivialLabels.size();
        for (int labelIndex = 0; labelIndex < labelCount; labelIndex++) {
            assignLayer(nonTrivialLabels.get(labelIndex), findPotentiallyWidestLayer(nonTrivialLabels, labelIndex));
        }
    }

    /**
     * Assigns all layers that either only have a single layer to be assigned to or that can be assigned to a layer
     * without increasing that layer's width. Returns the remaining labels in a list.
     */
    private List<LabelDummyInfo> performTrivialAssignments(final List<LabelDummyInfo> labelDummyInfos) {
        List<LabelDummyInfo> remainingLabels = new ArrayList<>(labelDummyInfos.size());
        
        for (LabelDummyInfo labelDummyInfo : labelDummyInfos) {
            if (labelDummyInfo.leftmostLayerId == labelDummyInfo.rightmostLayerId) {
                // Assign to only available layer and remove from list
                assignLayer(labelDummyInfo, labelDummyInfo.leftmostLayerId);
                
            } else if (!assignToWiderLayer(labelDummyInfo)) {
                // Ending up here means that we didn't find a layer wide enough for the node
                remainingLabels.add(labelDummyInfo);
            }
        }
        
        return remainingLabels;
    }
    
    /**
     * Assigns the given label dummy to the first layer wide enough to house it. If there is no such layer, the label
     * remains unassigned.
     * 
     * @return {@code true} if we succeeded in assigning the label to a layer.
     */
    private boolean assignToWiderLayer(final LabelDummyInfo labelDummyInfo) {
        // Check if the label dummy can be assigned a layer that already is at least as wide as the dummy
        double dummyWidth = labelDummyInfo.labelDummy.getSize().x;
        List<Layer> validLayers = labelDummyInfo.labelDummy.getGraph().getLayers().subList(
                labelDummyInfo.leftmostLayerId,
                labelDummyInfo.rightmostLayerId + 1);
        
        for (Layer layer : validLayers) {
            if (layer.getSize().x >= dummyWidth) {
                assignLayer(labelDummyInfo, layer.id);
                return true;
            }
        }
        
        // Ending up here means that we didn't find a layer wide enough for our label
        return false;
    }
    
    /**
     * Returns the index of the layer with the largest potential width.
     * 
     * @param labelDummyInfos
     *            label dummy infos, sorted descendingly by label width.
     * @param labelIndex
     *            index of the label dummy info for which we're currently looking for a layer. All labels with smaller
     *            index have already been assigned to layers.
     * @return index of the layer with the maximum potential width.
     */
    private int findPotentiallyWidestLayer(final List<LabelDummyInfo> labelDummyInfos, final int labelIndex) {
        assert labelIndex >= 0 && labelIndex < labelDummyInfos.size();
        
        int labelCount = labelDummyInfos.size();
        LabelDummyInfo labelDummyInfo = labelDummyInfos.get(labelIndex);
        double labelDummyWidth = labelDummyInfo.labelDummy.getSize().x;
        
        // Iterate over the label's valid layers
        int widestLayerIndex = labelDummyInfo.leftmostLayerId;
        double widestLayerWidth = 0;
        
        for (int layer = labelDummyInfo.leftmostLayerId; layer <= labelDummyInfo.rightmostLayerId; layer++) {
            // If the layer is already at least as large as the current label, simply return that layer
            if (labelDummyWidth <= layerWidths[layer]) {
                return layer;
            }
            
            // The initial potential width is less wide than the label (otherwise we would have already returned)
            double potentialWidth = layerWidths[layer];
            
            // Find the largest unassigned label that is part of this layer
            LabelDummyInfo largestUnassignedLabel = null;
            for (int label = labelIndex + 1; label < labelCount; label++) {
                // Check if the label can be placed in the current layer
                LabelDummyInfo currLabelInfo = labelDummyInfos.get(label);
                if (currLabelInfo.leftmostLayerId <= layer && currLabelInfo.rightmostLayerId >= layer) {
                    largestUnassignedLabel = currLabelInfo;
                }
            }
            
            // Update layer's potential size
            if (largestUnassignedLabel != null) {
                potentialWidth = Math.max(potentialWidth, largestUnassignedLabel.labelDummy.getSize().x);
            }
            
            // Update widest layer (if there are multiple widest layers, we use the leftmost one)
            if (potentialWidth > widestLayerWidth) {
                widestLayerIndex = layer;
                widestLayerWidth = potentialWidth;
            }
        }
        
        return widestLayerIndex;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Swapping Utilities
    
    /**
     * Assigns the given label dummy to the layer with the given index. Updates the layer size if the label enlarges
     * the layer.
     */
    private void assignLayer(final LabelDummyInfo labelDummyInfo, final int targetLayerIndex) {
        assert targetLayerIndex < labelDummyInfo.labelDummy.getGraph().getLayers().size();
        
        // If the label dummy is not in the target layer yet, swap it with the long edge dummy that is
        if (targetLayerIndex != labelDummyInfo.leftmostLayerId + labelDummyInfo.leftLongEdgeDummies.size()) {
            swapNodes(labelDummyInfo.labelDummy,
                    labelDummyInfo.ithDummyNode(targetLayerIndex - labelDummyInfo.leftmostLayerId));
        }
        
        // Update the size information of the label dummy's new layer
        int newLayerId = labelDummyInfo.labelDummy.getLayer().id;
        layerWidths[newLayerId] = Math.max(layerWidths[newLayerId], labelDummyInfo.labelDummy.getSize().x);
        
        for (LLabel label : labelDummyInfo.labelDummy.getProperty(InternalProperties.REPRESENTED_LABELS)) {
            label.setProperty(INCLUDE_LABEL, true);
        }
    }

    /**
     * Swaps the two given label dummy with the given long edge dummy. The dummies are assumed to be part of the same
     * long edge. The label dummy's new layer's width in the {@link #layerWidths} array is enlarged if the label dummy
     * is wider than the layer currently is.
     */
    private void swapNodes(final LNode labelDummy, final LNode longEdgeDummy) {
        // Find the layers and the positions inside the layers of the dummy nodes. We need the positions later since
        // we run after crossing minimization and have to keep the order of nodes the same. An alternative for this
        // method would be not to change the layers and connections of the two nodes but to switch all of their
        // properties instead, but we reckon that might actually be more work
        Layer layer1 = labelDummy.getLayer();
        Layer layer2 = longEdgeDummy.getLayer();
        
        int dummy1LayerPosition = layer1.getNodes().indexOf(labelDummy);
        int dummy2LayerPosition = layer2.getNodes().indexOf(longEdgeDummy);
        
        // Detect incoming and outgoing ports of the nodes (this of course assumes that there's just one of each kind,
        // which should be true for long edge and label dummy nodes)
        LPort inputPort1 = labelDummy.getPorts(PortType.INPUT).iterator().next();
        LPort outputPort1 = labelDummy.getPorts(PortType.OUTPUT).iterator().next();
        LPort inputPort2 = longEdgeDummy.getPorts(PortType.INPUT).iterator().next();
        LPort outputPort2 = longEdgeDummy.getPorts(PortType.OUTPUT).iterator().next();
        
        // Store incoming and outgoing edges
        LEdge[] incomingEdges1 = LGraphUtil.toEdgeArray(inputPort1.getIncomingEdges());
        LEdge[] outgoingEdges1 = LGraphUtil.toEdgeArray(outputPort1.getOutgoingEdges());
        LEdge[] incomingEdges2 = LGraphUtil.toEdgeArray(inputPort2.getIncomingEdges());
        LEdge[] outgoingEdges2 = LGraphUtil.toEdgeArray(outputPort2.getOutgoingEdges());

        // Put first dummy into second dummy's layer and reroute second dummy's edges to first dummy
        labelDummy.setLayer(dummy2LayerPosition, layer2);
        for (LEdge edge : incomingEdges2) {
            edge.setTarget(inputPort1);
        }
        for (LEdge edge : outgoingEdges2) {
            edge.setSource(outputPort1);
        }

        // Put second dummy into first dummy's layer and reroute first dummy's edges to second dummy
        longEdgeDummy.setLayer(dummy1LayerPosition, layer1);
        for (LEdge edge : incomingEdges1) {
            edge.setTarget(inputPort2);
        }
        for (LEdge edge : outgoingEdges1) {
            edge.setSource(outputPort2);
        }
    }
    
    /**
     * Updates the {@link InternalProperties#LONG_EDGE_BEFORE_LABEL_DUMMY} property of long edge dummy
     * nodes preceding and succeeding the given label dummy node.
     */
    private void updateLongEdgeSourceLabelDummyInfo(final LabelDummyInfo labelDummyInfo) {
        // Predecessors
        doUpdateLongEdgeLabelDummyInfo(
                labelDummyInfo.labelDummy,
                node -> node.getIncomingEdges().iterator().next().getSource().getNode(),
                true);
        
        // We may want to do things to the successors as well at some point
    }
    
    /**
     * Does the actual work of setting the long edge source or target of all nodes before or after
     * the given label dummy node to {@code null}.
     * 
     * @param labelDummy
     *            the label dummy node to start from.
     * @param nextElement
     *            a function that, given a node, returns the node to process next. Use this to
     *            decide whether to target all successors or all predecessors of the label dummy
     *            node.
     * @param value
     *            the new property value.
     */
    private void doUpdateLongEdgeLabelDummyInfo(final LNode labelDummy,
            final Function<LNode, LNode> nextElement, final boolean value) {
        
        LNode longEdgeDummy = nextElement.apply(labelDummy);
        while (longEdgeDummy.getType() == NodeType.LONG_EDGE) {
            longEdgeDummy.setProperty(InternalProperties.LONG_EDGE_BEFORE_LABEL_DUMMY, value);
            longEdgeDummy = nextElement.apply(longEdgeDummy);
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Label Dummy Info Class
    
    /**
     * Encapsulates information about a label dummy along with the long edge dummies to its left and to its right.
     */
    private static final class LabelDummyInfo {
        
        /** The label dummy node. */
        private LNode labelDummy;
        /** The label placement strategy to be used with this label. */
        private CenterEdgeLabelPlacementStrategy placementStrategy = null;
        /** The long edge dummies to the left of the label dummy. May well be empty. */
        private List<LNode> leftLongEdgeDummies = new ArrayList<>();
        /** The long edge dummies to the right of the label dummy. May well be empty. */
        private List<LNode> rightLongEdgeDummies = new ArrayList<>();
        /** ID of the leftmost layer the label dummy's long edge has a dummy node in. */
        private int leftmostLayerId;
        /** ID of the rightmost layer the label dummy's long edge has a dummy node in. */
        private int rightmostLayerId;
        
        /**
         * Create a new instance to represent the given label dummy. Automatically collects long edge dummies to its
         * left and to its right.
         * 
         * @param labelDummy
         *            the label dummy.
         */
        private LabelDummyInfo(final LNode labelDummy,
                final CenterEdgeLabelPlacementStrategy defaultPlacementStrategy) {
            
            this.labelDummy = labelDummy;
            this.placementStrategy = defaultPlacementStrategy;
            
            // Gather long edge dummies that are part of the label dummy's edge and find the end layer IDs
            gatherLeftLongEdgeDummies();
            gatherRightLongEdgeDummies();
            
            if (leftLongEdgeDummies.isEmpty()) {
                leftmostLayerId = labelDummy.getLayer().id;
            } else {
                leftmostLayerId = leftLongEdgeDummies.get(0).getLayer().id;
            }
            
            if (rightLongEdgeDummies.isEmpty()) {
                rightmostLayerId = labelDummy.getLayer().id;
            } else {
                rightmostLayerId = rightLongEdgeDummies.get(rightLongEdgeDummies.size() - 1).getLayer().id;
            }

            // Check if the label wants to deviate from the default placement strategy
            for (LLabel label : labelDummy.getProperty(InternalProperties.REPRESENTED_LABELS)) {
                // Take the first override we can find
                if (label.hasProperty(LayeredOptions.EDGE_LABELS_CENTER_LABEL_PLACEMENT_STRATEGY)) {
                    placementStrategy = label.getProperty(LayeredOptions.EDGE_LABELS_CENTER_LABEL_PLACEMENT_STRATEGY);
                    break;
                }
            }
        }

        /**
         * Gathers the long edge dummies left of the label dummy which are part of its long edge. The order they are
         * in reflects the order of layers, that is, the label dummy is the successor of the last long edge dummy in
         * the list.
         */
        private void gatherLeftLongEdgeDummies() {
            LNode source = labelDummy;
            do {
                source = source.getIncomingEdges().iterator().next().getSource().getNode();
                if (source.getType() == NodeType.LONG_EDGE) {
                    leftLongEdgeDummies.add(source);
                }
            } while (source.getType() == NodeType.LONG_EDGE);
            
            // The list is currently not in the order we would expect, so produce a reversed version
            leftLongEdgeDummies = Lists.reverse(leftLongEdgeDummies);
        }
        
        /**
         * Gathers the long edge dummies right of the label dummy which are part of its long edge. The order they are
         * in reflects the order of layers, that is, the label dummy is the predecessor of the first long edge dummy in
         * the list.
         */
        private void gatherRightLongEdgeDummies() {
            LNode target = labelDummy;
            do {
                target = target.getOutgoingEdges().iterator().next().getTarget().getNode();
                if (target.getType() == NodeType.LONG_EDGE) {
                    rightLongEdgeDummies.add(target);
                }
            } while (target.getType() == NodeType.LONG_EDGE);
        }
        
        /**
         * Returns the total number of long edge and label dummy nodes along the label dummy's edge.
         */
        public int totalDummyCount() {
            return rightmostLayerId - leftmostLayerId + 1;
        }
        
        /**
         * Returns the i-th dummy node from the edge through the long edge dummies and the label dummy.
         */
        private LNode ithDummyNode(final int i) {
            if (i < leftLongEdgeDummies.size()) {
                // The i-th dummy is a long edge dummy to the label dummy's left
                return leftLongEdgeDummies.get(i);
                
            } else if (i == leftLongEdgeDummies.size()) {
                // The i-th dummy is the label dummy itself
                return labelDummy;
                
            } else {
                // the i-th dummy is a long edge dummy to the label dummy's right
                return rightLongEdgeDummies.get(i - leftLongEdgeDummies.size() - 1);
            }
        }
    }
    
}
