/*******************************************************************************
 * Copyright (c) 2009, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.PriorityQueue;
import java.util.Queue;

import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.BoxLayouterOptions;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.graph.ElkNode;

import com.google.common.collect.Lists;

/**
 * A layout algorithm that does not take edges into account, but treats all nodes as isolated boxes.
 * This is useful for parts of a diagram that consist of objects without connections, such as
 * parallel regions in Statecharts.
 * 
 * <p>
 * The box layouter allows to configure the underlying packing strategy, see the {@link PackingMode} enumeration.
 * </p>
 * 
 * <p>
 * TODO there's duplicate code in this class, handling either the simple mode or any of the grouped modes. 
 * It should be possible to re-use the grouped code for the simple packing mode though.
 * </p>
 * 
 * <p>
 * MIGRATE The code needs to be changed to apply padding to the coordinates immediately. 
 * </p>
 * 
 * <p>
 * MIGRATE The box layout provider does not support hyperedges yet.
 * </p>
 */
public class BoxLayoutProvider extends AbstractLayoutProvider {
    
    /** default value for aspect ratio. */
    public static final double DEF_ASPECT_RATIO = 1.3;

    @Override
    public void layout(final ElkNode layoutNode, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Box layout", 2);
        
        float objSpacing = layoutNode.getProperty(BoxLayouterOptions.SPACING_NODE_NODE).floatValue();
        ElkPadding padding = layoutNode.getProperty(BoxLayouterOptions.PADDING);
        boolean expandNodes = layoutNode.getProperty(BoxLayouterOptions.EXPAND_NODES);
        boolean interactive = layoutNode.getProperty(BoxLayouterOptions.INTERACTIVE);

        switch (layoutNode.getProperty(BoxLayouterOptions.BOX_PACKING_MODE)) {
        case SIMPLE:
            placeBoxes(layoutNode, objSpacing, padding, expandNodes, interactive);
            break;

        default:
            // any of the groups
            placeBoxesGrouping(layoutNode, objSpacing, padding, expandNodes);
        }

        progressMonitor.done();
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Simple Mode

    /**
     * Place the boxes in the given parent node. The children are first sorted appropriately. Furthermore, adjust the
     * parent size to fit the bounding box.
     * 
     * @param parentNode parent node
     * @param objSpacing minimal spacing between elements
     * @param borderSpacing spacing to the border
     * @param expandNodes if true, the nodes are expanded to fill their parent
     * @param interactive whether position should be considered instead of size
     */
    private void placeBoxes(final ElkNode parentNode, final double objSpacing, final ElkPadding padding,
            final boolean expandNodes, final boolean interactive) {
        
        List<ElkNode> sortedBoxes = sort(parentNode, interactive);

        // Work on a copy of the minimum size to avoid changing the property's value
        KVector minSize = ElkUtil.effectiveMinSizeConstraintFor(parentNode);
        
        Double aspectRatio = parentNode.getProperty(BoxLayouterOptions.ASPECT_RATIO);
        if (aspectRatio == null || aspectRatio <= 0) {
            aspectRatio = DEF_ASPECT_RATIO;
        }

        // do place the boxes
        KVector parentSize = placeBoxes(sortedBoxes, objSpacing, padding, minSize.x, minSize.y,
                expandNodes, aspectRatio);

        // adjust parent size
        ElkUtil.resizeNode(parentNode, parentSize.x, parentSize.y, false, true);
    }

    /**
     * Sorts nodes according to priority and size or position. Nodes with higher priority are
     * put to the start of the list.
     * 
     * @param parentNode parent node
     * @param interactive whether position should be considered instead of size
     * @return sorted list of children
     */
    private List<ElkNode> sort(final ElkNode parentNode, final boolean interactive) {
        List<ElkNode> sortedBoxes = new ArrayList<>(parentNode.getChildren());
        
        Collections.sort(sortedBoxes, new Comparator<ElkNode>() {
            public int compare(final ElkNode child1, final ElkNode child2) {
                Integer prio1 = child1.getProperty(BoxLayouterOptions.PRIORITY);
                if (prio1 == null) {
                    prio1 = 0;
                }
                
                Integer prio2 = child2.getProperty(BoxLayouterOptions.PRIORITY);
                if (prio2 == null) {
                    prio2 = 0;
                }
                
                if (prio1 > prio2) {
                    return -1;
                } else if (prio1 < prio2) {
                    return 1;
                } else {
                    // boxes have same priority - compare their position or size
                    if (interactive) {
                        int c = Double.compare(child1.getY(), child2.getY());
                        if (c != 0) {
                            return c;
                        }
                        c = Double.compare(child1.getX(), child2.getX());
                        if (c != 0) {
                            return c;
                        }
                    }
                    double size1 = child1.getWidth() * child1.getHeight();
                    double size2 = child2.getWidth() * child2.getHeight();
                    return Double.compare(size1, size2);
                }
            }
        });
        
        return sortedBoxes;
    }

    /**
     * Place the boxes of the given sorted list according to their order in the list.
     * 
     * @param sortedBoxes sorted list of boxes
     * @param minSpacing minimal spacing between elements
     * @param borderSpacing spacing to the border
     * @param minTotalWidth minimal width of the parent node
     * @param minTotalHeight minimal height of the parent node
     * @param expandNodes if true, the nodes are expanded to fill their parent
     * @param aspectRatio the desired aspect ratio
     * @return the bounding box of the resulting layout
     */
    private KVector placeBoxes(final List<ElkNode> sortedBoxes, final double minSpacing,
            final ElkPadding padding, final double minTotalWidth, final double minTotalHeight,
            final boolean expandNodes, final double aspectRatio) {
        
        // determine the maximal row width by the maximal box width and the total area
        double maxRowWidth = 0.0f;
        double totalArea = 0.0f;
        for (ElkNode box : sortedBoxes) {
            ElkUtil.resizeNode(box);
            maxRowWidth = Math.max(maxRowWidth, box.getWidth());
            totalArea += box.getWidth() * box.getHeight();
        }
        
        // calculate std deviation
        //  rationale: the greater the diversity of box sizes, the more space will be 'wasted',
        //  contributing to the total area. We address this by adding x*n*stddev to the area.
        //  TODO x should be assessed empirically, for the moment set it to 1
        double mean = totalArea / sortedBoxes.size();
        double stddev = areaStdDev(sortedBoxes, mean);
        
        totalArea += (sortedBoxes.size() * 1 * stddev);
        // add top and bottom padding to total area
        totalArea += Math.sqrt(totalArea) * (padding.getBottom() + padding.getTop());
        // add left and right padding
        totalArea += Math.sqrt(totalArea) * (padding.getRight());


        // calculate the required row width w to achieve the desired aspect ratio,
        //  i.e.:  w*h=area s.t. w/h=dar  ->  w=sqrt(area * dar) 
        maxRowWidth = Math.max(maxRowWidth, Math.sqrt(totalArea * aspectRatio)) + padding.getLeft();

        // place nodes iteratively into rows
        double xpos = padding.getLeft();
        double ypos = padding.getTop();
        double highestBox = 0.0f;
        double broadestRow = padding.getHorizontal();
        LinkedList<Integer> rowIndices = new LinkedList<>();
        rowIndices.add(Integer.valueOf(0));
        LinkedList<Double> rowHeights = new LinkedList<>();
        ListIterator<ElkNode> boxIter = sortedBoxes.listIterator();
        while (boxIter.hasNext()) {
            ElkNode box = boxIter.next();
            double width = box.getWidth();
            double height = box.getHeight();
            if (xpos + width > maxRowWidth) {
                // place box into the next row
                if (expandNodes) {
                    rowHeights.addLast(Double.valueOf(highestBox));
                    rowIndices.addLast(Integer.valueOf(boxIter.previousIndex()));
                }
                xpos = padding.getLeft();
                ypos += highestBox + minSpacing;
                highestBox = 0.0f;
                broadestRow = Math.max(broadestRow, padding.getHorizontal() + width);
            }
            box.setLocation(xpos, ypos);
            broadestRow = Math.max(broadestRow, xpos + width + padding.getRight());
            highestBox = Math.max(highestBox, height);
            xpos += width + minSpacing;
        }
        broadestRow = Math.max(broadestRow, minTotalWidth);
        double totalHeight = ypos + highestBox + padding.getBottom();
        if (totalHeight < minTotalHeight) {
            highestBox += minTotalHeight - totalHeight;
            totalHeight = minTotalHeight;
        }

        // expand nodes if required
        if (expandNodes) {
            xpos = padding.getLeft();
            boxIter = sortedBoxes.listIterator();
            rowIndices.addLast(Integer.valueOf(sortedBoxes.size()));
            ListIterator<Integer> rowIndexIter = rowIndices.listIterator();
            int nextRowIndex = rowIndexIter.next();
            rowHeights.addLast(Double.valueOf(highestBox));
            ListIterator<Double> rowHeightIter = rowHeights.listIterator();
            double rowHeight = 0.0f;
            while (boxIter.hasNext()) {
                if (boxIter.nextIndex() == nextRowIndex) {
                    xpos = padding.getLeft();
                    rowHeight = rowHeightIter.next();
                    nextRowIndex = rowIndexIter.next();
                }
                ElkNode box = boxIter.next();
                double oldHeight = box.getHeight();
                box.setHeight(rowHeight);
                double newHeight = rowHeight;
                if (boxIter.nextIndex() == nextRowIndex) {
                    double newWidth = broadestRow - xpos - padding.getRight();
                    double oldWidth = box.getWidth();
                    box.setWidth(newWidth);
                    // Content alignment
                    ElkUtil.translate(box, new KVector(newWidth, newHeight), new KVector(oldWidth, oldHeight));
                }
                xpos += box.getWidth() + minSpacing;
            }
        }

        // return parent size
        return new KVector(broadestRow, totalHeight);
    }

    private double areaStdDev(final List<ElkNode> boxes, final double mean) {
        double variance = 0;
        for (ElkNode box : boxes) {
            variance += Math.pow(box.getWidth() * box.getHeight() - mean, 2);
        }
        double stddev = Math.sqrt(variance / (boxes.size() - 1));
        return stddev;
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Grouping Modes
    
    /**
     * Place the boxes of the given sorted list according to their order in the list.
     * Furthermore, adjust the parent size to fit the bounding box.
     * 
     * @param parentNode parent node
     * @param objSpacing minimal spacing between elements
     * @param borderSpacing spacing to the border
     * @param expandNodes if true, the nodes are expanded to fill their parent
     */
    private void placeBoxesGrouping(final ElkNode parentNode,
            final float objSpacing, final ElkPadding padding, final boolean expandNodes) {
        
        // Work on a copy of the minimum size to avoid changing the property's value
        KVector minSize = new KVector(parentNode.getProperty(BoxLayouterOptions.NODE_SIZE_MINIMUM));
        minSize.x = Math.max(minSize.x - padding.getLeft() - padding.getRight(), 0);
        minSize.y = Math.max(minSize.y - padding.getTop() - padding.getBottom(), 0);
        
        Double aspectRatio = parentNode.getProperty(BoxLayouterOptions.ASPECT_RATIO);
        if (aspectRatio == null || aspectRatio <= 0) {
            aspectRatio = DEF_ASPECT_RATIO;
        }

        // wrap boxes in groups
        List<Group> groups = Lists.newArrayList();
        for (ElkNode node : parentNode.getChildren()) {
            Group g = new Group(node);
            groups.add(g);
        }

        // pack according to the mode (sorting is done within every method)
        PackingMode mode = parentNode.getProperty(BoxLayouterOptions.BOX_PACKING_MODE);
        List<Group> toBePlaced;
        // no border spacing for the grouped nodes!
        switch (mode) {
        case GROUP_INC:
            toBePlaced = mergeAndPlaceInc(groups, objSpacing, minSize.x, minSize.y, expandNodes, aspectRatio);
            break;
        case GROUP_DEC:
            toBePlaced = mergeAndPlaceDec(groups, objSpacing, minSize.x, minSize.y, expandNodes, aspectRatio);
            break;
        default: // GROUP_MIXED
            toBePlaced = mergeAndPlaceMixed(groups, objSpacing, minSize.x, minSize.y, expandNodes, aspectRatio);
            break;
        }
        
        Group finalGroup = new Group(toBePlaced);
        // do final packing of the 'out-most' boxes
        KVector parentSize = placeInnerBoxes(finalGroup, objSpacing, padding, minSize.x, minSize.y,
                expandNodes, aspectRatio);

        // adjust parent size
        ElkUtil.resizeNode(parentNode, parentSize.x, parentSize.y, false, true);
    }
    
    /**
     * Place the boxes of the given sorted list according to their order in the list.
     * 
     * @param sortedBoxes sorted list of boxes
     * @param minSpacing minimal spacing between elements
     * @param borderSpacing spacing to the border
     * @param minTotalWidth minimal width of the parent node
     * @param minTotalHeight minimal height of the parent node
     * @param expandNodes if true, the nodes are expanded to fill their parent
     * @param aspectRatio the desired aspect ratio
     * @return the bounding box of the resulting layout
     */
    private KVector placeInnerBoxes(final Group group, final double minSpacing,
            final ElkPadding padding, final double minTotalWidth, final double minTotalHeight,
            final boolean expandNodes, final double aspectRatio) {
        // determine the maximal row width by the maximal box width and the total area
        double maxRowWidth = 0.0f;
        double totalArea = 0.0f;
        for (Group box : group.groups) {
            if (box.node != null) {
                ElkUtil.resizeNode(box.node);
            }
            maxRowWidth = Math.max(maxRowWidth, box.getWidth());
            totalArea += box.getWidth() * box.getHeight();
        }
        
        // calculate std deviation
        //  rationale: the greater the diversity of box sizes, the more space will be 'wasted',
        //  contributing to the total area. We address this by adding x*n*stddev to the area.
        //  TODO x should be assessed empirically, for the moment set it to 1 (GROUP_DEC worked better with 2 though)
        double mean = totalArea / group.groups.size();
        double stddev = areaStdDev2(group.groups, mean);
        
        double sdInfluence = 1;
        totalArea += (group.groups.size() * sdInfluence * stddev);

        // calculate the required row width w to achieve the desired aspect ratio,
        //  i.e.:  w*h=area s.t. w/h=dar  ->  w=sqrt(area * dar) 
        maxRowWidth = (double) Math.max(maxRowWidth, Math.sqrt(totalArea * aspectRatio)) + padding.getLeft();
        
        // place nodes iteratively into rows
        double xpos = padding.getLeft();
        double ypos = padding.getTop();
        double highestBox = 0.0f;
        double broadestRow = padding.getHorizontal();
        LinkedList<Integer> rowIndices = new LinkedList<Integer>();
        rowIndices.add(Integer.valueOf(0));
        LinkedList<Double> rowHeights = new LinkedList<>();
        ListIterator<Group> boxIter = group.groups.listIterator();
        Group last = null;
        List<Group> bottoms = Lists.newArrayList();
        while (boxIter.hasNext()) {
            Group box = boxIter.next();
            double width = box.getWidth();
            double height = box.getHeight();
            if (xpos + width > maxRowWidth) {
                // place box into the next row
                if (expandNodes) {
                    rowHeights.addLast(Double.valueOf(highestBox));
                    rowIndices.addLast(Integer.valueOf(boxIter.previousIndex()));
                    group.right.add(last);
                    bottoms.clear();
                }
                xpos = padding.getLeft();
                ypos += highestBox + minSpacing;
                highestBox = 0.0f;
                broadestRow = Math.max(broadestRow, padding.getHorizontal() + width);
            }
            bottoms.add(box);
            box.translate(xpos, ypos);
            broadestRow = Math.max(broadestRow, xpos + width + padding.getRight());
            highestBox = Math.max(highestBox, height);
            xpos += width + minSpacing;
            last = box;
        }
        group.bottom.addAll(bottoms);
        group.right.add(bottoms.get(bottoms.size() - 1));
        broadestRow = Math.max(broadestRow, minTotalWidth);
        double totalHeight = ypos + highestBox + padding.getBottom();
        if (totalHeight < minTotalHeight) {
            highestBox += minTotalHeight - totalHeight;
            totalHeight = minTotalHeight;
        }

        // expand nodes if required
        if (expandNodes) {
            xpos = padding.getLeft();
            boxIter = group.groups.listIterator();
            rowIndices.addLast(Integer.valueOf(group.groups.size()));
            ListIterator<Integer> rowIndexIter = rowIndices.listIterator();
            int nextRowIndex = rowIndexIter.next();
            rowHeights.addLast(Double.valueOf(highestBox));
            ListIterator<Double> rowHeightIter = rowHeights.listIterator();
            double rowHeight = 0.0f;
            while (boxIter.hasNext()) {
                if (boxIter.nextIndex() == nextRowIndex) {
                    xpos = padding.getLeft();
                    rowHeight = rowHeightIter.next();
                    nextRowIndex = rowIndexIter.next();
                }
                Group box = boxIter.next();
                box.setHeight(rowHeight);
                if (boxIter.nextIndex() == nextRowIndex) {
                    double newWidth = broadestRow - xpos - padding.getRight();
                    double oldWidth = box.getWidth();
                    box.setWidth(newWidth);
                    box.translateInnerNodes((newWidth - oldWidth) / 2, 0.0f);
                }
                xpos += box.getWidth() + minSpacing;
            }
        }

        // return parent size
        return new KVector(broadestRow, totalHeight);
    }
    
    private double areaStdDev2(final List<Group> boxes, final double mean) {
        double variance = 0;
        for (Group box : boxes) {
            variance += Math.pow(box.area() - mean, 2);
        }
        double stddev = (double) Math.sqrt(variance / (boxes.size() - 1));
        return stddev;
    }

    // Now the packing modes ///////////////////////////////////////////////////////////////////////////

    /**
     * Start with large nodes and collect smaller nodes afterwards until the large node's 
     * half size is reached, then groups the small nodes.
     */
    private List<Group> mergeAndPlaceDec(final List<Group> groups, 
            final double objSpacing,
            final double minWidth, final double minHeight, 
            final boolean expandNodes,
            final double aspectRatio) {
        
        // sort in decreasing area
        // the code below makes relies on this
        Collections.sort(groups, (g1, g2) -> -Double.compare(g1.area(), g2.area()));
        
        Queue<Group> boxQueue = Lists.newLinkedList(groups);
        List<Group> toBePlaced = Lists.newArrayList();
        List<Group> maybeGroup = Lists.newArrayList();
        
        // current 'large' node and the combined area of 'smaller' nodes
        Group boxToBeat = null;         
        double collectedArea = 0;
        
        while (!boxQueue.isEmpty()) {
            Group box = boxQueue.poll();
            if ((boxToBeat == null)
                    // prevent 'box' from getting its own 'innerGroup'
                    || (boxToBeat.area() / 2 < box.area())) {
                boxToBeat = box;
                toBePlaced.add(box);
            } else {
                collectedArea += box.area();
                maybeGroup.add(box);
                
                if (maybeGroup.size() > 1 
                        && ((collectedArea > boxToBeat.area() / 2)
                        // make sure the remaining nodes in the queue get grouped
                        || boxQueue.isEmpty())) {
                    
                    // combine
                    Group innerGroup = new Group(maybeGroup);
                    // TODO there may be a better choice for the inner aspect ratio
                    double innerAspectRatio = boxToBeat.getWidth() / boxToBeat.getHeight();
                    KVector groupSize = placeInnerBoxes(innerGroup, objSpacing, new ElkPadding(), minWidth, minHeight,
                            expandNodes, innerAspectRatio);
                    innerGroup.size.reset().add(groupSize);
                    
                    boxToBeat = innerGroup;
                    toBePlaced.add(innerGroup);

                    collectedArea = 0;
                    maybeGroup.clear();
                } 
            }
        }
        
        // add remaining boxes
        toBePlaced.addAll(maybeGroup);
        return toBePlaced;
    }
    
    /**
     * Use a priority queue, once a large node is encountered, 
     * check if there's an (arbitrary) number of small nodes
     * that can be combined to come close to the larger node's size. 
     */
    private List<Group> mergeAndPlaceMixed(final List<Group> groups, 
            final double objSpacing,
            final double minWidth, final double minHeight, 
            final boolean expandNodes,
            final double aspectRatio) {

        double[] cumAreaArray = new double[groups.size()];
        
        PriorityQueue<Group> pq =
                new PriorityQueue<BoxLayoutProvider.Group>((g1, g2) -> Double.compare(g1.area(), g2.area()));
        pq.addAll(groups);
        
        int index = 0;
        List<Group> toBePlaced = Lists.newArrayList();
        
        while (!pq.isEmpty()) {
            // only take a look, we might not want to remove it
            Group box = pq.peek();
            
            // we need at least two boxes in the 'toBePlaced' list that can be merged 
            // thus 'index > 1'
            if (index > 1 && (box.area() / 2f > cumAreaArray[0])) {
                
                // TODO this should be a binary search 
                int anIndex = 0;
                // find the largest index
                // we can use at most 'toBePlaced.size()' groups 
                while (anIndex < toBePlaced.size() - 1 && box.area() / 2f > cumAreaArray[anIndex]) {
                    anIndex++;
                }
                
                // place the selected groups
                List<Group> select = toBePlaced.subList(0, anIndex + 1);
                Group innerGroup = new Group(select);
                // TODO there's probably a better strategy to select this aspect ratio
                double innerAspectRatio = box.getWidth() / box.getHeight();

                KVector groupSize = placeInnerBoxes(innerGroup, objSpacing, new ElkPadding(), minWidth, minHeight,
                        expandNodes, innerAspectRatio);
                innerGroup.size.reset().add(groupSize);
                
                // now start over
                pq.add(innerGroup);
                List<Group> remain = toBePlaced.subList(anIndex + 1, toBePlaced.size());
                pq.addAll(remain);
                toBePlaced.clear();
                index = 0;
                Arrays.fill(cumAreaArray, 0);
                
            } else {
                // ok, now remove the group from the queue
                pq.poll();
                if (index > 0) {
                    cumAreaArray[index] = cumAreaArray[index - 1];
                }
                cumAreaArray[index] += box.area();
                index++;
                toBePlaced.add(box);
            }
        }
        
        return toBePlaced;
    }
    
    /**
     * In increasing node area order: pick nodes until the next node's half size is larger than
     * the current cumulated area. Problematic for certain graphs, where 
     * the cumulation becomes quite large early. 
     */
    private List<Group> mergeAndPlaceInc(final List<Group> groups, 
            final double objSpacing, 
            final double minWidth, final double minHeight, 
            final boolean expandNodes,
            final double aspectRatio) {

        // sort increasingly 
        Collections.sort(groups, (g1, g2) -> Double.compare(g1.area(), g2.area()));
        
        ListIterator<Group> groupIterator = groups.listIterator();
        List<Group> toBePlaced = Lists.newArrayList();
        double commonArea = 0;

        while (groupIterator.hasNext()) {
            Group g = groupIterator.next();
            if (!toBePlaced.isEmpty() && g.area() > (commonArea * 2)) {
                // merge the current set of groups into a common group
                Group merged = new Group(toBePlaced);
                double innerAspectRatio = g.getWidth() / g.getHeight();
                
                // place the nodes in the group
                KVector groupSize = placeInnerBoxes(merged, objSpacing, new ElkPadding(), minWidth, minHeight,
                        expandNodes, innerAspectRatio);
                merged.size.reset().add(groupSize);

                // reset
                toBePlaced.clear();
                commonArea = 0;

                // add the merged group and the current group to
                // the active list 
                toBePlaced.add(merged);
                toBePlaced.add(g);
                commonArea = merged.area() + g.area();
            } else {
                // add the group to the active set
                toBePlaced.add(g);
                commonArea += g.area();
            }
        }
        
        return toBePlaced;
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Classes and Enumerations
    
    /**
     * Configures the packing mode used by the {@link BoxLayoutProvider}.
     * If {@value #SIMPLE} is not required, {@link #GROUP_DEC} could improve the packing and decrease the area. 
     * {@value #GROUP_MIXED} and {@link #GROUP_INC} may, in very specific scenarios, work better.
     */
    public enum PackingMode {
        
        /** In order to use either of {@link BoxLayouterOptions#PRIORITY} 
         * {@link BoxLayouterOptions#INTERACTIVE}, this mode must be used.*/
        SIMPLE,
        
        /** Tries to group boxes iterating the boxes in decreasing area. */
        GROUP_DEC,

        /**
         * Iterating in increasing area, tries to be smart by grouping a subset of small nodes such that their combined
         * area is close to half of the area of a larger node.
         */
        GROUP_MIXED,
        
        /** Tries to group boxes iterating the boxes in increasing area. */
        GROUP_INC,       
    }
    
    /**
     * Internal representation of a group. A group either wraps a single node, or holds a collection of groups.
     * In essence it delegates position/dimension changes to the inner-most groups (wrapping nodes).
     */
    private static class Group {
        
        /** The node represented by this group, may be {@code null}. */
        private ElkNode node;
        /** The set of groups wrapped in this group, may only be used if {@code node} is {@code null}. */
        private List<Group> groups; 
        /** If this group represents other groups, {@code size} represents the bounding rectangle. */
        private KVector size;

        /** Stores the bottom-most nodes of this group, used if {@link CoreOptions#EXPAND_NODES} is active .*/
        private List<Group> bottom;
        /** Stores the right-most nodes of this group, used if {@link CoreOptions#EXPAND_NODES} is active .*/
        private List<Group> right;
        
        Group(final ElkNode node) {
            this.node = node;
            node.setLocation(0, 0);
        }
        
        Group(final Iterable<Group> groups) {
            this.groups = Lists.newArrayList(groups);
            this.bottom = Lists.newArrayList();
            this.right = Lists.newArrayList();
            this.size = new KVector();
        }

        public double area() {
            return getWidth() * getHeight();
        }
        
        public double getWidth() {
            if (node != null) {
                return node.getWidth();
            }
            return size.x;
        }
        
        public double getHeight() {
            if (node != null) {
                return node.getHeight();
            }
            return size.y;
        }

        /**
         * Sets (increases) the width of this group. For node-wrapping groups this directly alters the
         * node's with. For groups of groups this only makes sense in conjunction with
         * {@link CoreOptions#EXPAND_NODES}, effectively increasing all inner-most nodes (by recursion) of the
         * {@link #right} collection as to fit their parent group's width.
         */
        public void setWidth(final double w) {
            assert w > getWidth();
            if (node != null) {
                node.setWidth(w);
            } else {
                // increase width of the rightmost nodes
                double delta = w - getWidth();
                for (Group g : right) {
                    g.setWidth(g.getWidth() + delta);
                }
            }
        }
        
        /**
         * Sets (increases) the height of this group. See {@link #setWidth(float)} for further details.
         */
        public void setHeight(final double h) {
            assert h > getHeight();
            if (node != null) {
                node.setHeight(h);
            } else {
                // increase bottom most groups
                double delta = h - getHeight();
                for (Group g : bottom) {
                    g.setHeight(g.getHeight() + delta);
                }
            }
        }
     
        /**
         * For groups wrapping nodes, the passed {@code x} and {@code y} are added to the node's
         * location. For group-wrapping groups, the method is called recursively. 
         * 
         * As opposed to {@link #translateInnerNodes(double, double)} this method alters the position of a wrapped
         * node's location.
         */
        public void translate(final double x, final double y) {
            if (node != null) {
                node.setX(node.getX() + x);
                node.setY(node.getY() + y);
            } else {
                for (Group g : groups) {
                    g.translate(x, y);
                }
            }
        }

        /**
         * Translates the nested nodes of wrapped hierarchical nodes.
         * 
         * As opposed to {@link #translate(double, double)} this method translates the nested nodes of a hierarchical
         * node. See also {@link ElkUtil#translate(ElkNode, double, double)}.
         */
        public void translateInnerNodes(final double x, final double y) {
            if (node != null) {
                ElkUtil.translate(node, x, y);
            } else {
                for (Group g : groups) {
                    g.translateInnerNodes(x, y);
                }
            }
        }
        
        @Override
        public String toString() {
            if (node != null) {
                return node.toString();
            } else {
                return groups.toString();
            }
        }
    }
    
}
