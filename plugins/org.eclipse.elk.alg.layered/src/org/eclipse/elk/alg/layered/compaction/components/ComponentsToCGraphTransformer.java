/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.layered.compaction.components;

import static org.eclipse.elk.core.options.PortSide.SIDES_EAST;
import static org.eclipse.elk.core.options.PortSide.SIDES_EAST_SOUTH;
import static org.eclipse.elk.core.options.PortSide.SIDES_EAST_SOUTH_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_EAST_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_EAST;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_EAST_SOUTH_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_EAST_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_SOUTH;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_SOUTH_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_SOUTH;
import static org.eclipse.elk.core.options.PortSide.SIDES_SOUTH_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_WEST;

import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.eclipse.elk.alg.layered.compaction.oned.CGraph;
import org.eclipse.elk.alg.layered.compaction.oned.CGroup;
import org.eclipse.elk.alg.layered.compaction.oned.CNode;
import org.eclipse.elk.alg.layered.compaction.oned.ICGraphTransformer;
import org.eclipse.elk.alg.layered.compaction.oned.ISpacingsHandler;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.Pair;

import java.util.EnumMap;

/**
 * Transforms {@link IConnectedComponents} into a representing {@link CGraph} which can be passed to
 * a {@link org.eclipse.elk.alg.layered.compaction.oned.OneDimensionalCompactor
 * OneDimensionalCompactor} to obtain a more compact drawing.
 * 
 * @param <N>
 *            the type of the contained nodes
 * @param <E>
 *            the type of the contained edges
 */
public class ComponentsToCGraphTransformer<N, E> implements
        ICGraphTransformer<IConnectedComponents<N, E>> {

    /** The {@link CGroup} resulting from this transformation. */
    private CGraph cGraph;
    
    private double spacing;

    // internal mappings 
    private Map<IComponent<N, E>, KVector> oldPosition = new HashMap<>();
    private Map<IComponent<N, E>, CRectNode> offsets = new HashMap<>();
    
    // a global offset and the new graph size after layout
    private KVector globalOffset;
    private KVector graphSize;

    // maps the external extension to the CNodes by which they are represented
    private Map<IExternalExtension<E>, Pair<CGroup, CNode>> externalExtensions = new HashMap<>();
    private Map<Direction, List<Pair<CGroup, CNode>>> externalPlaceholder = new EnumMap<>(Direction.class);
    
    /**
     * @param spacing required spacing between connected IComponents.
     */
    public ComponentsToCGraphTransformer(final double spacing) {
        this.spacing = spacing;
    }

    /**
     * @param c
     *            a IComponent
     * @return the offset to c's original position after layout.
     */
    public KVector getOffset(final IComponent<N, E> c) {
        KVector cOffset = oldPosition.get(c).clone().sub(offsets.get(c).rect.getPosition());
        return cOffset;
    }
    
    /**
     * @return the globalOffset
     */
    public KVector getGlobalOffset() {
        return globalOffset;
    }
    
    /**
     * @return the graphSize
     */
    public KVector getGraphSize() {
        return graphSize;
    }
    
    /**
     * @return the externalExtensions
     */
    public Map<IExternalExtension<E>, Pair<CGroup, CNode>> getExternalExtensions() {
        return externalExtensions;
    }
    
    /**
     * @return the externalPlaceholder
     */
    public Map<Direction, List<Pair<CGroup, CNode>>> getExternalPlaceholder() {
        return externalPlaceholder;
    }
    
    /* -----------------------------------------------------------
     *                    Graph Transformation
     * ----------------------------------------------------------- */
    
    @Override
    public CGraph transform(final IConnectedComponents<N, E> ccs) {

        cGraph = new CGraph(EnumSet.allOf(Direction.class));

        for (IComponent<N, E> comp : ccs.getComponents()) {

            CGroup group = new CGroup();
            cGraph.cGroups.add(group);

            // convert the hull of the graph's elements without external edges
            for (ElkRectangle rect : comp.getHull()) {

                CRectNode rectNode = new CRectNode(rect);
                setLock(rectNode, comp.getExternalExtensionSides());
                
                if (!oldPosition.containsKey(comp)) {
                    oldPosition.put(comp, new KVector(rect.x, rect.y));
                    offsets.put(comp, rectNode);
                }

                cGraph.cNodes.add(rectNode);
                group.addCNode(rectNode);
            }

            // now prepare the rectangles for the external extensions
            // they can be added to the CGraph on demand later on
            for (IExternalExtension<E> ee : comp.getExternalExtensions()) {
                CRectNode rectNode = new CRectNode(ee.getRepresentor());
                externalExtensions.put(ee, Pair.of(group, rectNode));
                setLock(rectNode, comp.getExternalExtensionSides());
                
                // placeholder
                if (ee.getPlaceholder() != null) {
                    CRectNode rectPlaceholder = new CRectNode(ee.getPlaceholder(), 1d);
                    setLock(rectPlaceholder, comp.getExternalExtensionSides());
                    CGroup dummyGroup = new CGroup();
                    dummyGroup.addCNode(rectPlaceholder);
                    externalPlaceholder.computeIfAbsent(ee.getDirection(), k -> new ArrayList<>())
                            .add(Pair.of(group, rectPlaceholder));
                }
            }
        }

        return cGraph;
    }
    
    /**
     * The locks should assure that a component c that is connected to external ports on the sides
     * WEST and SOUTH stays in the bottom left corner. Otherwise edges may become unnecessarily
     * long.
     */
    private void setLock(final CNode cNode, final Set<PortSide> portSides) {
        
        if (portSides.isEmpty()) {
            cNode.lock.set(true, true, true, true);            
        }
        
        // regulars
        if (portSides.equals(SIDES_NORTH)) {
            cNode.lock.set(true, true, true, false);
        }
        if (portSides.equals(SIDES_EAST)) {
            cNode.lock.set(false, true, true, true);
        }
        if (portSides.equals(SIDES_SOUTH)) {
            cNode.lock.set(true, true, false, true);
        }
        if (portSides.equals(SIDES_WEST)) {
            cNode.lock.set(true, false, true, true);
        }

        // corners
        if (portSides.equals(SIDES_NORTH_EAST)) {
            cNode.lock.set(false, true, true, false);
        }
        if (portSides.equals(SIDES_EAST_SOUTH)) {
            cNode.lock.set(false, true, false, true);
        }
        if (portSides.equals(SIDES_SOUTH_WEST)) {
            cNode.lock.set(true, false, false, true);
        }
        if (portSides.equals(SIDES_NORTH_WEST)) {
            cNode.lock.set(true, false, true, false);
        }

        // opposite
        if (portSides.equals(SIDES_NORTH_SOUTH)) {
            cNode.lock.set(true, true, true, true);
        }
        if (portSides.equals(SIDES_EAST_WEST)) {
            cNode.lock.set(true, true, true, true);
        }

        // triumvirates
        if (portSides.equals(SIDES_NORTH_SOUTH)) {
            cNode.lock.set(true, true, true, true);
        }
        if (portSides.equals(SIDES_EAST_SOUTH_WEST)) {
            cNode.lock.set(true, true, true, true);
        }
        if (portSides.equals(SIDES_NORTH_SOUTH_WEST)) {
            cNode.lock.set(true, true, true, true);
        }
        if (portSides.equals(SIDES_NORTH_EAST_WEST)) {
            cNode.lock.set(true, true, true, true);
        }

        // all of them!
        if (portSides.equals(SIDES_NORTH_EAST_SOUTH_WEST)) {
            cNode.lock.set(true, true, true, true);
        }
    }

    /* -----------------------------------------------------------
     *                    Layout Application 
     * ----------------------------------------------------------- */
    
    @Override
    public void applyLayout() {

        for (CNode n : cGraph.cNodes) {
            n.applyElementPosition();
        }

        // calculating new graph size and offset
        KVector topLeft = new KVector(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        KVector bottomRight = new KVector(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        for (CNode cNode : cGraph.cNodes) {
            topLeft.x = Math.min(topLeft.x, cNode.hitbox.x);
            topLeft.y = Math.min(topLeft.y, cNode.hitbox.y);
            bottomRight.x = Math.max(bottomRight.x, cNode.hitbox.x + cNode.hitbox.width);
            bottomRight.y = Math.max(bottomRight.y, cNode.hitbox.y + cNode.hitbox.height);
        }
        
        // note that placeholdes "can move" during compaction and are _not_ fixed to the boundary
        // which allows to use them to determine the maximum position
        for (List<Pair<CGroup, CNode>> placeholders : externalPlaceholder.values()) {
            for (Pair<CGroup, CNode> placeholder : placeholders) {
                CNode cNode = placeholder.getSecond();
                topLeft.x = Math.min(topLeft.x, cNode.hitbox.x);
                topLeft.y = Math.min(topLeft.y, cNode.hitbox.y);
                bottomRight.x = Math.max(bottomRight.x, cNode.hitbox.x + cNode.hitbox.width);
                bottomRight.y = Math.max(bottomRight.y, cNode.hitbox.y + cNode.hitbox.height);
            }
        }
        
        globalOffset = topLeft.clone().negate();
        graphSize = bottomRight.clone().sub(topLeft);

        // resetting lists
        cGraph.cGroups.clear();
        cGraph.cNodes.clear();
    }

    /**
     * {@link CNode} implementation for {@link ElkRectangle}s representing 
     * parts of the rectangular convex hull of a connected IComponent.
     */
    private final class CRectNode extends CNode {

        private ElkRectangle rect;
        private Double individualSpacing;

        private CRectNode(final ElkRectangle rect) {
            this(rect, null);
        }
        
        private CRectNode(final ElkRectangle rect, final Double spacing) {
            this.rect = rect;
            this.hitbox = new ElkRectangle(rect.x, rect.y, rect.width, rect.height);
            this.individualSpacing = spacing;
        }

        @Override
        public double getHorizontalSpacing() {
            return individualSpacing != null ? individualSpacing : spacing;
        }
        
        @Override
        public double getVerticalSpacing() {
           return individualSpacing != null ? individualSpacing : spacing;
        }
        
        @Override
        public void applyElementPosition() {
            rect.x = hitbox.x;
            rect.y = hitbox.y;
        }

        @Override
        public double getElementPosition() {
            return rect.x;
        }

        @Override
        public String toString() {
            // for debug output be quiet! :)
            return "";
            // return "CRectNode " + hitbox;
        }
    }
    
    /**
     * Component compaction-specific spacing handler that can cope with placeholders.
     * During compaction a placeholder can be moved right next to an adjacent rectangle;
     * no spacing has to be preserved.
     */
    public static final ISpacingsHandler<? super CNode> SPACING_HANDLER = new ISpacingsHandler<CNode>() {
        @Override
        public double getHorizontalSpacing(final CNode cNode1, final CNode cNode2) {
            return Math.min(cNode1.getHorizontalSpacing(), cNode2.getHorizontalSpacing());
        }

        @Override
        public double getVerticalSpacing(final CNode cNode1, final CNode cNode2) {
            return Math.min(cNode1.getVerticalSpacing(), cNode2.getVerticalSpacing());
        }
    };
}
