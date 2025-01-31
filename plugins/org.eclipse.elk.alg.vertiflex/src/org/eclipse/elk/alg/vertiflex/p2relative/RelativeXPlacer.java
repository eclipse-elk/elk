/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.vertiflex.p2relative;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.eclipse.elk.alg.vertiflex.InternalProperties;
import org.eclipse.elk.alg.vertiflex.VertiFlexLayoutPhases;
import org.eclipse.elk.alg.vertiflex.options.VertiFlexOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.graph.ElkNode;

/**
 * Node placer that positions nodes horizontally using coordinates that are relative to parent nodes.
 *
 */
public class RelativeXPlacer implements ILayoutPhase<VertiFlexLayoutPhases, ElkNode> {

    private double spacingNodeNode;
    private boolean considerNodeModelOrder;

    // a constant for moving every Outline to a minimal y-pos
    private static final double MINIMAL_Y = -100.0;
    
    @Override
    public void process(final ElkNode graph, final IElkProgressMonitor progressMonitor) {

        progressMonitor.begin("XPlacer", 1);
        
        spacingNodeNode = graph.getProperty(CoreOptions.SPACING_NODE_NODE);
        considerNodeModelOrder = graph.getProperty(VertiFlexOptions.CONSIDER_NODE_MODEL_ORDER);
        
        if (!graph.getChildren().isEmpty()) {
            ElkNode parent = graph.getProperty(InternalProperties.ROOT_NODE);
            
            switch (graph.getProperty(VertiFlexOptions.LAYOUT_STRATEGY)) {
                case STRAIGHT:
                    recursiveStraightlinePlacement(parent);
                    break;
                case BEND:
                    recursiveBentlinePlacement(parent);
                    break;
                default:
                    break;   
            }
        }

        progressMonitor.done();
    }
    
    /** Computes the distance between two outlines. */
    private double outlineDistance(final OutlineNode outline1, final OutlineNode outline2) {
        
        OutlineNode changedOutline1 = new OutlineNode(outline1.getRelativeX(), MINIMAL_Y, 
                new OutlineNode(0.0, outline1.getAbsoluteY(), outline1.getNext()));
        OutlineNode changedOutline2 = new OutlineNode(outline2.getRelativeX(), MINIMAL_Y, 
                new OutlineNode(0.0, outline2.getAbsoluteY(), outline2.getNext()));
        
        // the return value
        double dist = changedOutline1.getRelativeX() - changedOutline2.getRelativeX();
        
        OutlineNode o1;
        OutlineNode o2;
        double x1;
        double x2;
        double deltaX, deltaY;
        double newdist;
        
        
        // first run (compare points of o1  with o2)
        o1 = changedOutline1;
        o2 = changedOutline2;
        x1 = o1.getRelativeX();
        x2 = o2.getRelativeX();
        while (o1 != null && !o2.isLast()) {
            if (o2.getNext().getAbsoluteY() > o1.getAbsoluteY()) {
             // now we compare
                deltaX = o2.getNext().getRelativeX();
                deltaY = o2.getNext().getAbsoluteY() - o2.getAbsoluteY();
                newdist = x1 - x2 - ((o1.getAbsoluteY() - o2.getAbsoluteY()) * deltaX) / deltaY;
                
                dist = Math.max(dist, newdist);
                
                // now change o1
                o1 = o1.getNext();
                if (o1 != null) {
                    x1 += o1.getRelativeX();
                }
            } else {
                o2 = o2.getNext();
                x2 += o2.getRelativeX();
            }
        }
        
        // second run (compare points of o2  with o1)
        o1 = changedOutline1;
        o2 = changedOutline2;
        x1 = o1.getRelativeX();
        x2 = o2.getRelativeX();
        while (o2 != null && !o1.isLast()) {
            if (o1.getNext().getAbsoluteY() > o2.getAbsoluteY()) {
                // now we compare
                deltaX = o1.getNext().getRelativeX();
                deltaY = o1.getNext().getAbsoluteY() - o1.getAbsoluteY();
                newdist = x1 - x2 + ((o2.getAbsoluteY() - o1.getAbsoluteY()) * deltaX) / deltaY;
                
                dist = Math.max(dist, newdist);
                
                // now change o2
                o2 = o2.getNext();
                if (o2 != null) {
                    x2 += o2.getRelativeX();
                }
            } else {
                o1 = o1.getNext();
                x1 += o1.getRelativeX();
            }
      
        }
        
        return dist;
    }
    
    /**
     * This is the recursive function that calculates the layout for one node and it's children.
     * Children are placed such that straight edges can later be drawn from their parent to each of them.
     * 
     * If ConsiderModelOrder is set to false, all children are arranged in a semi-circle with the parent
     * initially positioned above the lowest child and then shifted toward the center of the children as 
     * far as possible.
     * 
     * If ConsiderModelOrder is set to true, groups of children are arranged such that they maintain their
     * inherent model order. Furthermore, the reading direction is oriented left to right and top to bottom 
     * as far as possible without violating the vertical position constraints and the straight edge routing
     * requirement.
     */
    private void recursiveStraightlinePlacement(final ElkNode graph) {
        
        makeSimpleOutlines(graph);
        
        if (!graph.getOutgoingEdges().isEmpty()) {
            
            // get all children
            List<ElkNode> children = new ArrayList<>();
            for (int i = 0; i < graph.getOutgoingEdges().size(); i++) {
                ElkNode child = (ElkNode) graph.getOutgoingEdges().get(i).getTargets().get(0);
                recursiveStraightlinePlacement(child);
                children.add(child);
            }
            // now the children of this node get sorted to form a semi-circle. This allows routing straight edges
            // without overlaps while keeping the layout relatively compact.
            sortSubTrees(children);
            
            // now the children get stuffed together, using the outlines.
            for (int i = 0; i < children.size() - 1; i++) {
                bundleChildren(children.get(0), children.get(i), children.get(i + 1));
            }
            
            // now we need to move the root to the middle of the nodes.
            // we calculate the point of the child with the lowest y-position to avoid overlapping.
            // if there is more than one lowest child, the root will be positioned in the middle of them.
                int pos = 0;
                double maxDepth = 0.0;
                int maxDepthStartPos = 0;
                while (pos < children.size() && children.get(pos).getY() >= maxDepth) {
                    if (children.get(pos).getY() > maxDepth) {
                        maxDepthStartPos = pos;
                        maxDepth = children.get(pos).getY();
                    }
                    pos += 1;
                }
                double moveRoot = 0.0;
                if (pos > 0) {
                    moveRoot = (children.get(maxDepthStartPos).getX() + children.get(pos - 1).getX()) / 2.0 
                            - graph.getX();
                }
            
            if (!graph.getProperty(VertiFlexOptions.CONSIDER_NODE_MODEL_ORDER)) {
                double betterMoveRoot = (children.get(0).getX() + children.get(children.size() - 1).getX() 
                        + children.get(children.size() - 1).getWidth() - graph.getWidth()) / 2.0 - graph.getX();
                double newMoveRoot;
    
                if (betterMoveRoot < moveRoot) {
                    OutlineNode rightOutline;
                    double rightOutlineX, posX;
                    for (int i = 0; i < maxDepthStartPos; i++) {
                        for (int j = i + 1; j < maxDepthStartPos + 1; j++) {
                            rightOutline = children.get(i).getProperty(InternalProperties.RIGHT_OUTLINE);
                            rightOutlineX = children.get(i).getX() + rightOutline.getRelativeX();
                            posX = children.get(j).getX() + children.get(j).getWidth() / 2.0;
                            while (rightOutline != null && rightOutline.getAbsoluteY() < maxDepth) {
    
                                newMoveRoot = posX - graph.getWidth() / 2.0 + (posX - rightOutlineX) * ((graph.getY() 
                                        + graph.getHeight()) - maxDepth) / (maxDepth - rightOutline.getAbsoluteY());
                                betterMoveRoot = Math.max(betterMoveRoot, newMoveRoot);
    
                                rightOutline = rightOutline.getNext();
                                if (rightOutline != null) {
                                    rightOutlineX += rightOutline.getRelativeX();
                                }
                            }
                        }
                    }
                    moveRoot = betterMoveRoot;
                }
    
                if (betterMoveRoot > moveRoot) {
                    OutlineNode leftOutline;
                    double leftOutlineX, posX;
                    for (int i = pos; i < children.size(); i++) {
                        for (int j = pos - 1; j < i; j++) {
                            leftOutline = children.get(i).getProperty(InternalProperties.LEFT_OUTLINE);
                            leftOutlineX = children.get(i).getX() + leftOutline.getRelativeX();
                            posX = children.get(j).getX() + children.get(j).getWidth() / 2.0;
                            while (leftOutline != null && leftOutline.getAbsoluteY() < maxDepth) {
    
                                newMoveRoot = posX - graph.getWidth() / 2.0 + (posX - leftOutlineX) * ((graph.getY() 
                                        + graph.getHeight()) - maxDepth) / (maxDepth - leftOutline.getAbsoluteY());
                                betterMoveRoot = Math.min(betterMoveRoot, newMoveRoot);
    
                                leftOutline = leftOutline.getNext();
                                if (leftOutline != null) {
                                    leftOutlineX += leftOutline.getRelativeX();
                                }
                            }
                        }
                    }
                    moveRoot = betterMoveRoot;
                }
            }
            
            for (ElkNode child: children) {
                child.setX(child.getX() - moveRoot);
            }
            
            
            double newX;
            // left outline update
            OutlineNode graphLeftOutline = graph.getProperty(InternalProperties.LEFT_OUTLINE);
            OutlineNode leftChildOutline = children.get(0).getProperty(InternalProperties.LEFT_OUTLINE);
            
            newX = children.get(0).getX() + leftChildOutline.getRelativeX()  - graphLeftOutline.getRelativeX();
            graphLeftOutline.getNext().getNext().getNext().setNext(
                    new OutlineNode(newX, leftChildOutline.getAbsoluteY(), leftChildOutline.getNext()));
            
            // right outline update
            OutlineNode graphRightOutline = graph.getProperty(InternalProperties.RIGHT_OUTLINE);
            OutlineNode rightChildOutline = children.get(children.size() - 1)
                    .getProperty(InternalProperties.RIGHT_OUTLINE);
            
            newX = children.get(children.size() - 1).getX() 
                    + rightChildOutline.getRelativeX() - graphRightOutline.getRelativeX();
            graphRightOutline.getNext().getNext().getNext().setNext(
                    new OutlineNode(newX, rightChildOutline.getAbsoluteY(), rightChildOutline.getNext()));
            
            // update outlineMaxY
            // update min und max for x and y
            for (ElkNode child: children) {
                graph.setProperty(InternalProperties.OUTLINE_MAX_DEPTH, Math.max(graph.getProperty(InternalProperties.
                        OUTLINE_MAX_DEPTH), child.getProperty(InternalProperties.OUTLINE_MAX_DEPTH)));
                graph.setProperty(InternalProperties.MIN_X, Math.min(graph.getProperty(InternalProperties.MIN_X), 
                        child.getX() + child.getProperty(InternalProperties.MIN_X)));
                graph.setProperty(InternalProperties.MAX_X, Math.max(graph.getProperty(InternalProperties.MAX_X), 
                        child.getX() + child.getProperty(InternalProperties.MAX_X)));
            }
            graph.setProperty(InternalProperties.MAX_Y, graph.getProperty(InternalProperties.OUTLINE_MAX_DEPTH));
            
        }
    }

    /** 
     * Place nodes while maintaining model order and computing bendpoints for the later edges. 
     * The model order is fully kept intact and bendpoints for the edges are computed  so that 
     * overlap free edge routing is still possible.
     */
    private void recursiveBentlinePlacement(final ElkNode graph) {
        
        // set up initial outlines for all nodes
        makeSimpleOutlines(graph);
        
        // termination condition for the recursion, while a node has children continue
        if (!graph.getOutgoingEdges().isEmpty()) {
            
            // get all children
            List<ElkNode> children = new ArrayList<>();
            for (int i = 0; i < graph.getOutgoingEdges().size(); i++) {
                ElkNode child = (ElkNode) graph.getOutgoingEdges().get(i).getTargets().get(0);
                recursiveBentlinePlacement(child);
                children.add(child);
            }
            
            int childrenSize = children.size();
            
            
            // now the children get stuffed together, using the outlines.
            for (int i = 0; i < children.size() - 1; i++) {
                bundleChildren(children.get(0), children.get(i), children.get(i + 1));
            }
            
            double moveRoot = (children.get(0).getX() + children.get(0).getWidth() / 2.0 
                    + children.get(children.size() - 1).getX() + children.get(children.size() - 1).getWidth() / 2.0 
                    - graph.getWidth()) / 2.0 - graph.getX();
            
            for (ElkNode child: children) {
                child.setX(child.getX() - moveRoot);
                child.setProperty(InternalProperties.EDGE_BEND_HEIGHT, 
                        child.getProperty(InternalProperties.LEFT_OUTLINE).getAbsoluteY());
            }
            
            // set bendHeights for children right of the parent
            int i = 0;
            while (i < childrenSize - 1 && children.get(i).getX() + children.get(i).getWidth() 
                    + children.get(i).getProperty(CoreOptions.MARGINS).right - graph.getWidth() / 2.0 <= 0.0) {
                i++;
            }
            
            double globalBendHeight = children.get(i).getProperty(InternalProperties.EDGE_BEND_HEIGHT);
            for (int a = 0; a < childrenSize; a++) {
                if (globalBendHeight < children.get(a).getProperty(InternalProperties.EDGE_BEND_HEIGHT)) {
                    children.get(a).setProperty(InternalProperties.EDGE_BEND_HEIGHT, globalBendHeight);
                } else {
                    globalBendHeight = children.get(a).getProperty(InternalProperties.EDGE_BEND_HEIGHT);
                }
            }
            
            // set bendHeights for children left of the parent
            i = childrenSize - 1;
            while (i > 0 && children.get(i).getX() - children.get(i).getProperty(CoreOptions.MARGINS).left 
                    - graph.getWidth() / 2.0 >= 0.0) {
                i--;
            }
            
            if (i < childrenSize) {
                for (int a = i; a >= 0; a--) {
                    
                    if (globalBendHeight < children.get(a).getProperty(InternalProperties.EDGE_BEND_HEIGHT)) {
                        children.get(a).setProperty(InternalProperties.EDGE_BEND_HEIGHT, globalBendHeight);
                    } else {
                        globalBendHeight = children.get(a).getProperty(InternalProperties.EDGE_BEND_HEIGHT);
                    }
                }
            }
            
            double newX;
            OutlineNode newOutlinepart;
            // left outline update
            OutlineNode leftChildOutline = children.get(0).getProperty(InternalProperties.LEFT_OUTLINE);
            OutlineNode graphLeftOutline = graph.getProperty(InternalProperties.LEFT_OUTLINE);
            
            newX = children.get(0).getX() + leftChildOutline.getRelativeX() 
                    - graphLeftOutline.getRelativeX();
            newOutlinepart = new OutlineNode(0.0, leftChildOutline
                    .getAbsoluteY(), leftChildOutline.getNext());
            graphLeftOutline.getNext().getNext().getNext()
                    .setNext(new OutlineNode(newX, children.get(0).getProperty(InternalProperties.EDGE_BEND_HEIGHT), 
                            newOutlinepart));
            
            // right outline update
            OutlineNode graphRightOutline = graph.getProperty(InternalProperties.RIGHT_OUTLINE);
            OutlineNode rightChildOutline = children.get(childrenSize - 1).getProperty(InternalProperties.RIGHT_OUTLINE);
            
            newX = children.get(childrenSize - 1).getX() + rightChildOutline.getRelativeX() - graphRightOutline.getRelativeX();
            newOutlinepart = new OutlineNode(0.0, rightChildOutline.getAbsoluteY(), rightChildOutline.getNext());
            graphRightOutline.getNext().getNext().getNext().setNext(
                    new OutlineNode(
                            newX, children.get(childrenSize - 1).getProperty(InternalProperties.EDGE_BEND_HEIGHT), newOutlinepart
                    ));
            
            // update outlineMaxY
            // update min und max for x and y
            for (ElkNode child: children) {
                graph.setProperty(InternalProperties.OUTLINE_MAX_DEPTH, Math.max(graph.getProperty(InternalProperties
                        .OUTLINE_MAX_DEPTH), child.getProperty(InternalProperties.OUTLINE_MAX_DEPTH)));
                graph.setProperty(InternalProperties.MIN_X, Math.min(graph.getProperty(InternalProperties.MIN_X), 
                        child.getX() + child.getProperty(InternalProperties.MIN_X)));
                graph.setProperty(InternalProperties.MAX_X, Math.max(graph.getProperty(InternalProperties.MAX_X), 
                        child.getX() + child.getProperty(InternalProperties.MAX_X)));
            }
            graph.setProperty(InternalProperties.MAX_Y, graph.getProperty(InternalProperties.OUTLINE_MAX_DEPTH));
            
        }
        
    }
    
    /** Sorts the subTrees in a semi-circle. */
    private void sortSubTrees(final List<ElkNode> children) {

        // first, we sort the SubTrees by the Y-coordinate of their root.
        Collections.sort(children, new NodeComparator(false));

        List<ElkNode> a = new ArrayList<>();
        List<ElkNode> b = new ArrayList<>();

        if (considerNodeModelOrder) {
            splitNodesWithModelOrder(children, a, b);
        } else {
            // now we need to put them in a V-shape 
            // the deepest element gets ignored in the calculation of the widths.
            a.add(children.get(children.size() - 1));
            double widthA = 0.0, widthB = 0.0;
            for (int i = 1; i < children.size(); i++) {
                if (widthA <= widthB) {
                    a.add(children.get(children.size() - 1 - i));
                    widthA += children.get(children.size() - 1 - i).getWidth();
                } else {
                    b.add(children.get(children.size() - 1 - i));
                    widthB += children.get(children.size() - 1 - i).getWidth();
                }
            }
            Collections.reverse(a);
        }

        Collections.sort(b, new NodeComparator(true));
        a.addAll(b);
        for (int i = 0; i < children.size(); i++) {
            children.set(i, a.get(i));
        }

    }

    /** Split nodes into two lists, while maintaining a sensible model order for nodes on the same height. */
    private void splitNodesWithModelOrder(List<ElkNode> original, List<ElkNode> left, List<ElkNode> right) {
        // identify next subgroup (all at same y)
        //split group in it in the weighted middle (always add first half to left group, and second half to right group)
        if (original.size() == 0) {
            return;
        }
        if (original.size() == 1) {
            left.add(original.get(0));
        }
        if (original.size() == 2) {
            // use model order to decide which node to put into which list regardless what the height ordering says
            ElkNode first = original.get(0);
            ElkNode second = original.get(1);
            
            if (first.getProperty(InternalProperties.NODE_MODEL_ORDER) 
                    > second.getProperty(InternalProperties.NODE_MODEL_ORDER)) {
                left.add(second);
                right.add(first);
            } else {
                left.add(first);
                right.add(second);
            }
        }

        List<ElkNode> currentGroup = new ArrayList<>();
        Pair<Double, Double> widthLeftRight = new Pair<>();
        widthLeftRight.setFirst(0.0);
        widthLeftRight.setSecond(0.0);
        ElkNode current = original.get(0);
        currentGroup.add(current);
        for (int i = 1; i < original.size(); i++) {
            ElkNode next = original.get(i);
            if (Double.compare(current.getY(), next.getY()) == 0) {
                // add next node to group because it has the same height as the previous element
                currentGroup.add(next);
            } else {
                // split and add entire group
                int finalIndexOfLeft = splitGroup(currentGroup, widthLeftRight);
                left.addAll(currentGroup.subList(0, finalIndexOfLeft + 1));
                right.addAll(currentGroup.subList(finalIndexOfLeft + 1, currentGroup.size()));
                
                // reset group
                currentGroup = new ArrayList<>();
                currentGroup.add(next);
            }
            current = next;
            
            // if next element is last element split and add the current group now
            if (i == original.size() - 1) {
                int finalIndexOfLeft = splitGroup(currentGroup, widthLeftRight);
                left.addAll(currentGroup.subList(0, finalIndexOfLeft + 1));
                right.addAll(currentGroup.subList(finalIndexOfLeft + 1, currentGroup.size()));
            }
        }
    }
    
    /** Find a sensible splitting point to divide elements into a left and right list according to a given width of both
     * lists.
     * @param group the group to be split
     * @param widthLeftRight the current widths of the left and right lists
     * @return the index position at which to split the group
     */
    private int splitGroup(List<ElkNode> group, Pair<Double, Double> widthLeftRight) {
        if (group.size() == 1) {
            return 0;
        } else {
            double widthLeft = widthLeftRight.getFirst();
            double widthRight = widthLeftRight.getSecond();
            
            double totalNewWidth = 0;
            for (ElkNode node : group) {
                totalNewWidth += node.getWidth();
            }
            
            // linear equations to find desired widths to fill both lists so that they reach equal length
            //  I: a + x = b + y
            // II: x + y = t
            // where a and b are widthLeft and widthRight, x and y are the desired additions to left and right
            // and t is the totalNewWidth
            // solution: x = (t-a+b)/2, y = t - (t-a+b)/2
            
            double desiredLeft = (totalNewWidth - widthLeft + widthRight) / 2;
//            double desiredRight = totalNewWidth - desiredLeft;
            
            // add nodes to left side until desiredLeft is exceeded
            // remove last element and compare both solutions, take the one which is closer to the desired solution
            int i = 0;
            double newLeftWidth = 0;
            while (desiredLeft < newLeftWidth && i < group.size()) {
                newLeftWidth += group.get(i).getWidth();
            }
            
            double exceed = newLeftWidth - desiredLeft;
            double under = desiredLeft - (newLeftWidth - group.get(i).getWidth());
            
            int resultIndex;
            if (exceed > under) {
                resultIndex = i;
            } else {
                resultIndex = i - 1;
            }
            
            // add new widths to accumulated widths
            widthLeftRight.setFirst(widthLeft + exceed);
            double newRightWidth = 0;
            for (int j = resultIndex + 1; j < group.size(); j++) {
                newRightWidth += group.get(j).getWidth();
            }
            widthLeftRight.setSecond(widthRight + newRightWidth);
            
            return resultIndex;
        }
    }
    
    /** Create the initial outlines around a node. */
    private void makeSimpleOutlines(final ElkNode graph) {
        ElkMargin margins = graph.getProperty(CoreOptions.MARGINS);
        
        // set the properties for left and right outlines
        OutlineNode endpart;
        endpart = new OutlineNode(0.0, graph.getY() + graph.getHeight() + margins.bottom + spacingNodeNode / 2, 
                new OutlineNode(graph.getWidth() / 2.0, 
                        graph.getY() + graph.getHeight() + margins.bottom + spacingNodeNode / 2, null));
        
        graph.setProperty(InternalProperties.LEFT_OUTLINE, new OutlineNode((-margins.left - spacingNodeNode / 2) 
                + graph.getWidth() / 2.0, graph.getY() - margins.top - spacingNodeNode / 2, 
                new OutlineNode(-graph.getWidth() / 2.0, graph.getY() - margins.top, endpart)));
        
        endpart = new OutlineNode(0.0, graph.getY() + graph.getHeight() + margins.bottom,
                new OutlineNode(-graph.getWidth() / 2.0, 
                        graph.getY() + graph.getHeight() + margins.bottom + spacingNodeNode / 2, null));
        
        graph.setProperty(InternalProperties.RIGHT_OUTLINE, new OutlineNode(graph.getWidth() / 2.0 
                + margins.right + spacingNodeNode / 2, graph.getY() - margins.top,
                new OutlineNode(graph.getWidth() / 2.0, graph.getY() - margins.top - spacingNodeNode / 2, endpart)));
        
        // set min and max values
        graph.setProperty(InternalProperties.MIN_X, graph.getX() - margins.left);
        graph.setProperty(InternalProperties.MAX_X, 
                graph.getX() + margins.right + graph.getWidth());
        graph.setProperty(InternalProperties.MIN_Y, graph.getY() - margins.top);
        graph.setProperty(InternalProperties.MAX_Y,
                graph.getY() + margins.bottom + graph.getHeight());
        graph.setProperty(InternalProperties.OUTLINE_MAX_DEPTH, graph.getProperty(InternalProperties.LEFT_OUTLINE)
                .getNext().getNext().getAbsoluteY());
        
    }
    
    /**
     * Combines the individual outlines of several sibling nodes to form the outline of the tree that they belong to.
     */
    private void bundleChildren(final ElkNode leftSubtree, final ElkNode a, final ElkNode b) {
        
        double deltaX, deltaY, change;
        
        // calculate distance between the two parts
        double dist = outlineDistance(a.getProperty(InternalProperties.RIGHT_OUTLINE), 
                b.getProperty(InternalProperties.LEFT_OUTLINE));
        b.setX(a.getX() + dist);
        
        // enhance the left outline
        if (leftSubtree.getProperty(InternalProperties.OUTLINE_MAX_DEPTH) 
                < b.getProperty(InternalProperties.OUTLINE_MAX_DEPTH)) {
            OutlineNode lastL = leftSubtree.getProperty(InternalProperties.LEFT_OUTLINE);
            double lAbsX = lastL.getRelativeX() + leftSubtree.getX();
            
            // move to the end of leftSubtree
            while (!lastL.isLast()) {
                lastL = lastL.getNext();
                lAbsX += lastL.getRelativeX();
            }
            // find fitting position in the left outline of b
            OutlineNode bItterator = new OutlineNode(b.getProperty(InternalProperties.LEFT_OUTLINE).getRelativeX(),
                    MINIMAL_Y, b.getProperty(InternalProperties.LEFT_OUTLINE).getNext());
            double rAbsX = bItterator.getRelativeX() + b.getX();
            while (bItterator.getNext().getAbsoluteY() <= lastL.getAbsoluteY()) {
                bItterator = bItterator.getNext();
                rAbsX += bItterator.getRelativeX();
            }
            // now we calculate the change
            deltaX = bItterator.getNext().getRelativeX();
            deltaY = bItterator.getNext().getAbsoluteY() - bItterator.getAbsoluteY();
            change = ((lastL.getAbsoluteY() - bItterator.getAbsoluteY()) * deltaX) / deltaY;
            
            // now we calculate the new points
            double newX = -lAbsX + rAbsX + change;
            OutlineNode newNext = new OutlineNode(bItterator.getNext().getRelativeX() - change, bItterator.getNext()
                    .getAbsoluteY(), bItterator.getNext().getNext());
            lastL.setNext(new OutlineNode(newX, lastL.getAbsoluteY(), newNext));
            
            // now update outline_max_depth
            leftSubtree.setProperty(InternalProperties.OUTLINE_MAX_DEPTH, 
                    b.getProperty(InternalProperties.OUTLINE_MAX_DEPTH));
        }
        
        // enhance the right outline
        if (b.getProperty(InternalProperties.OUTLINE_MAX_DEPTH) < a.getProperty(InternalProperties.OUTLINE_MAX_DEPTH)) {
            OutlineNode lastB = b.getProperty(InternalProperties.RIGHT_OUTLINE);
            double rAbsX = lastB.getRelativeX() + b.getX();
            
            // move to the end of b
            while (!lastB.isLast()) {
                lastB = lastB.getNext();
                rAbsX += lastB.getRelativeX();
            }
            // find fitting position in the right outline of a
            OutlineNode aItterator = new OutlineNode(a.getProperty(InternalProperties.RIGHT_OUTLINE).getRelativeX(), 
                    MINIMAL_Y, a.getProperty(InternalProperties.RIGHT_OUTLINE).getNext());
            double aAbsX = aItterator.getRelativeX() + a.getX();
            while (aItterator.getNext().getAbsoluteY() <= lastB.getAbsoluteY()) {
                aItterator = aItterator.getNext();
                aAbsX += aItterator.getRelativeX();
            }
            // now we calculate the change
            deltaX = aItterator.getNext().getRelativeX();
            deltaY = aItterator.getNext().getAbsoluteY() - aItterator.getAbsoluteY();
            change = ((lastB.getAbsoluteY() - aItterator.getAbsoluteY()) * deltaX) / deltaY;
            
            // now we calculate the new points
            double newX = aAbsX - rAbsX + change;
            OutlineNode newNext = new OutlineNode(aItterator.getNext().getRelativeX() - change, 
                    aItterator.getNext().getAbsoluteY(), aItterator.getNext().getNext());
            lastB.setNext(new OutlineNode(newX, lastB.getAbsoluteY(), newNext));
            // now update outline_max_depth
            b.setProperty(InternalProperties.OUTLINE_MAX_DEPTH, a.getProperty(InternalProperties.OUTLINE_MAX_DEPTH));
        }
    }

    @Override
    public LayoutProcessorConfiguration<VertiFlexLayoutPhases, ElkNode> 
    getLayoutProcessorConfiguration(final ElkNode graph) {
        return null;
    }

}
