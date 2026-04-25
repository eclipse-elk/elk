/*******************************************************************************
 * Copyright (c) 2012, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.components;

import static org.eclipse.elk.core.options.PortSide.SIDES_EAST;
import static org.eclipse.elk.core.options.PortSide.SIDES_EAST_SOUTH;
import static org.eclipse.elk.core.options.PortSide.SIDES_EAST_SOUTH_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_EAST_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_NONE;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_EAST;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_EAST_SOUTH;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_EAST_SOUTH_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_EAST_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_SOUTH;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_SOUTH_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_NORTH_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_SOUTH;
import static org.eclipse.elk.core.options.PortSide.SIDES_SOUTH_WEST;
import static org.eclipse.elk.core.options.PortSide.SIDES_WEST;

import java.util.Collection;
import java.util.List;
import java.util.ArrayList;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.EdgeRouting;


/**
 * A graph placer that tries to place the components of a graph with taking connections to external
 * ports into account. This graph placer should only be used if the constraints applying to the
 * external ports are either {@code FREE} or {@code FIXED_SIDES}.
 * 
 * <p>This placer first greedily builds a list of {@link ComponentGroup} instances. It is greedy in
 * that it places a component in the first group it finds that is able to hold it. Afterwards, the
 * components in each group are placed. The different groups are then placed along a diagonal from
 * the top-left to the bottom-right corner.</p>
 * 
 * <p>The target graph must not be contained in the list of components.</p>
 */
class ComponentGroupGraphPlacer extends AbstractGraphPlacer {
    
    ///////////////////////////////////////////////////////////////////////////////
    // Variables
    
    /**
     * List of component groups holding the different components.
     */
    protected final List<ComponentGroup> componentGroups = new ArrayList<>();
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // AbstractGraphPlacer

    public void combine(final List<LGraph> components, final LGraph target) {
        componentGroups.clear();
        assert !components.contains(target);
        target.getLayerlessNodes().clear();
        
        // Check if there are any components to be placed
        if (components.isEmpty()) {
            target.getSize().x = 0;
            target.getSize().y = 0;
            return;
        }
        
        // Set the graph properties
        LGraph firstComponent = components.get(0);
        target.copyProperties(firstComponent);
        
        // Construct component groups
        for (LGraph component : components) {
            addComponent(component);
        }
        
        // Place components in each group
        KVector offset = new KVector();
        double componentSpacing = firstComponent.getProperty(LayeredOptions.SPACING_COMPONENT_COMPONENT);
        
        for (ComponentGroup group : componentGroups) {
            // Place the components
            KVector groupSize = placeComponents(group, componentSpacing);
            offsetGraphs(group.getComponents(), offset.x, offset.y);
            
            // Compute the new offset
            offset.x += groupSize.x;
            offset.y += groupSize.y;
        }
        
        // Set the graph's new size (the component group sizes include additional spacing
        // on the right and bottom sides which we need to subtract at this point)
        target.getSize().x = offset.x - componentSpacing;
        target.getSize().y = offset.y - componentSpacing;
        
        // if compaction is desired, do so!cing;
        if (firstComponent.getProperty(LayeredOptions.COMPACTION_CONNECTED_COMPONENTS)
                // the compaction only supports orthogonally routed edges
                && firstComponent.getProperty(LayeredOptions.EDGE_ROUTING) == EdgeRouting.ORTHOGONAL) {

            // apply graph offsets (which we reset later on)
            // since the compaction works in a common coordinate system
            for (LGraph h : components) {
                offsetGraph(h, h.getOffset().x, h.getOffset().y);
            }

            ComponentsCompactor compactor = new ComponentsCompactor();
            compactor.compact(components, target.getSize(), componentSpacing);

            // the compaction algorithm places components absolutely,
            // therefore we have to use the final drawing's offset
            for (LGraph h : components) {
                h.getOffset().reset().add(compactor.getOffset());
            }

            // set the new (compacted) graph size
            target.getSize().reset().add(compactor.getGraphSize());
        }

        // finally move the components to the combined graph
        for (ComponentGroup group : componentGroups) {
            moveGraphs(target, group.getComponents(), 0, 0);
        }
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // Component Group Building
    
    /**
     * Adds the given component to the first component group that has place for it.
     * 
     * @param component the component to be placed.
     */
    protected void addComponent(final LGraph component) {
        // Check if one of the existing component groups has some place left
        for (ComponentGroup group : componentGroups) {
            if (group.add(component)) {
                return;
            }
        }
        
        // Create a new component group for the component
        componentGroups.add(new ComponentGroup(component));
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // Component Placement
    
    /* Component placement works as follows.
     * 
     * We first go through all the possible combinations of external port side connections. For each
     * combination, we compute a placement for the components with the given combination, remembering
     * the amount of space the placement takes up.
     * 
     * We then go through all the combinations again, this time moving the components so as not to
     * overlap other components.
     */
    
    
    /**
     * Computes a placement for the components in the given component group.
     * 
     * @param group the group whose components are to be placed.
     * @param spacing the amount of space to leave between two components.
     * @return the group's size.
     */
    protected KVector placeComponents(final ComponentGroup group, final double spacing) {
        
        // Determine the spacing between two components
        // Place the different sector components and remember the amount of space their placement uses.
        // In this phase, we pretend that no other components are in the component group.
        KVector sizeC = placeComponentsInRows(
                group.getComponents(SIDES_NONE), spacing);
        KVector sizeN = placeComponentsHorizontally(
                group.getComponents(SIDES_NORTH), spacing);
        KVector sizeS = placeComponentsHorizontally(
                group.getComponents(SIDES_SOUTH), spacing);
        KVector sizeW = placeComponentsVertically(
                group.getComponents(SIDES_WEST), spacing);
        KVector sizeE = placeComponentsVertically(
                group.getComponents(SIDES_EAST), spacing);
        KVector sizeNW = placeComponentsHorizontally(
                group.getComponents(SIDES_NORTH_WEST), spacing);
        KVector sizeNE = placeComponentsHorizontally(
                group.getComponents(SIDES_NORTH_EAST), spacing);
        KVector sizeSW = placeComponentsHorizontally(
                group.getComponents(SIDES_SOUTH_WEST), spacing);
        KVector sizeSE = placeComponentsHorizontally(
                group.getComponents(SIDES_EAST_SOUTH), spacing);
        KVector sizeWE = placeComponentsVertically(
                group.getComponents(SIDES_EAST_WEST), spacing);
        KVector sizeNS = placeComponentsHorizontally(
                group.getComponents(SIDES_NORTH_SOUTH), spacing);
        KVector sizeNWE = placeComponentsHorizontally(
                group.getComponents(SIDES_NORTH_EAST_WEST), spacing);
        KVector sizeSWE = placeComponentsHorizontally(
                group.getComponents(SIDES_EAST_SOUTH_WEST), spacing);
        KVector sizeWNS = placeComponentsVertically(
                group.getComponents(SIDES_NORTH_SOUTH_WEST), spacing);
        KVector sizeENS = placeComponentsVertically(
                group.getComponents(SIDES_NORTH_EAST_SOUTH), spacing);
        KVector sizeNESW = placeComponentsHorizontally(
                group.getComponents(SIDES_NORTH_EAST_SOUTH_WEST), spacing);
        
        // Find the maximum height of the three rows and the maximum width of the three columns the
        // component group is divided into (we're adding a fourth row for WE components and a fourth
        // column for NS components to make the placement easier later)
        double colLeftWidth = ElkMath.maxd(sizeNW.x, sizeW.x, sizeSW.x, sizeWNS.x);
        double colMidWidth = ElkMath.maxd(sizeN.x, sizeC.x, sizeS.x, sizeNESW.x);
        double colNsWidth = sizeNS.x;
        double colRightWidth = ElkMath.maxd(sizeNE.x, sizeE.x, sizeSE.x, sizeENS.x);
        double rowTopHeight = ElkMath.maxd(sizeNW.y, sizeN.y, sizeNE.y, sizeNWE.y);
        double rowMidHeight = ElkMath.maxd(sizeW.y, sizeC.y, sizeE.y, sizeNESW.y);
        double rowWeHeight = sizeWE.y;
        double rowBottomHeight = ElkMath.maxd(sizeSW.y, sizeS.y, sizeSE.y, sizeSWE.y);
        
        // With the individual placements computed, we now move the components to their final place,
        // taking the size of other component placements into account (the NW, NWE, and WNS components
        // stay at coordinates (0,0) and thus don't need to be moved around)
        offsetGraphs(group.getComponents(SIDES_NONE),
                colLeftWidth + colNsWidth,
                rowTopHeight + rowWeHeight);
        offsetGraphs(group.getComponents(SIDES_NORTH_EAST_SOUTH_WEST),
                colLeftWidth + colNsWidth,
                rowTopHeight + rowWeHeight);
        offsetGraphs(group.getComponents(SIDES_NORTH),
                colLeftWidth + colNsWidth,
                0.0);
        offsetGraphs(group.getComponents(SIDES_SOUTH),
                colLeftWidth + colNsWidth,
                rowTopHeight + rowWeHeight + rowMidHeight);
        offsetGraphs(group.getComponents(SIDES_WEST),
                0.0,
                rowTopHeight + rowWeHeight);
        offsetGraphs(group.getComponents(SIDES_EAST),
                colLeftWidth + colNsWidth + colMidWidth,
                rowTopHeight + rowWeHeight);
        offsetGraphs(group.getComponents(SIDES_NORTH_EAST),
                colLeftWidth + colNsWidth + colMidWidth,
                0.0);
        offsetGraphs(group.getComponents(SIDES_SOUTH_WEST),
                0.0,
                rowTopHeight + rowWeHeight + rowMidHeight);
        offsetGraphs(group.getComponents(SIDES_EAST_SOUTH),
                colLeftWidth + colNsWidth + colMidWidth,
                rowTopHeight + rowWeHeight + rowMidHeight);
        offsetGraphs(group.getComponents(SIDES_EAST_WEST),
                0.0,
                rowTopHeight);
        offsetGraphs(group.getComponents(SIDES_NORTH_SOUTH),
                colLeftWidth,
                0.0);
        offsetGraphs(group.getComponents(SIDES_EAST_SOUTH_WEST),
                0.0,
                rowTopHeight + rowWeHeight + rowMidHeight);
        offsetGraphs(group.getComponents(SIDES_NORTH_EAST_SOUTH),
                colLeftWidth + colNsWidth + colMidWidth,
                0.0);
        
        // Compute this component group's size
        KVector componentSize = new KVector();
        componentSize.x = ElkMath.maxd(
                colLeftWidth + colMidWidth + colNsWidth + colRightWidth,
                sizeWE.x,
                sizeNWE.x,
                sizeSWE.x);
        componentSize.y = ElkMath.maxd(
                rowTopHeight + rowMidHeight + rowWeHeight + rowBottomHeight,
                sizeNS.y,
                sizeWNS.y,
                sizeENS.y);
        
        return componentSize;
    }
    
    /**
     * Places the given collection of components along a horizontal line.
     * 
     * @param components the components to place.
     * @param spacing the amount of space to leave between two components.
     * @return the space used by the component placement, including spacing to the right and to the
     *         bottom of the components.
     */
    protected KVector placeComponentsHorizontally(final Collection<LGraph> components,
            final double spacing) {
        
        KVector size = new KVector();
        
        // Iterate over the components and place them
        for (LGraph component : components) {
            offsetGraph(component, size.x, 0.0);
            
            size.x += component.getSize().x + spacing;
            size.y = Math.max(size.y, component.getSize().y);
        }
        
        // Add vertical spacing, if necessary
        if (size.y > 0.0) {
            size.y += spacing;
        }
        
        return size;
    }
    
    /**
     * Places the given collection of components along a vertical line.
     * 
     * @param components the components to place.
     * @param spacing the amount of space to leave between two components.
     * @return the space used by the component placement.
     */
    protected KVector placeComponentsVertically(final Collection<LGraph> components,
            final double spacing) {
        
        KVector size = new KVector();
        
        // Iterate over the components and place them
        for (LGraph component : components) {
            offsetGraph(component, 0.0, size.y);
            
            size.y += component.getSize().y + spacing;
            size.x = Math.max(size.x, component.getSize().x);
        }
        
        // Add horizontal spacing, if necessary
        if (size.x > 0.0) {
            size.x += spacing;
        }
        
        return size;
    }
    
    /**
     * Place the given collection of components in multiple rows.
     * 
     * @param components the components to place.
     * @param spacing the amount of space to leave between two components.
     * @return the space used by the component placement.
     */
    protected KVector placeComponentsInRows(final Collection<LGraph> components,
            final double spacing) {
        
        /* This code is basically taken from the SimpleRowGraphPlacer. */
        
        // Check if there actually are any components
        if (components.isEmpty()) {
            return new KVector();
        }
        
        // Determine the maximal row width by the maximal box width and the total area
        double maxRowWidth = 0.0f;
        double totalArea = 0.0f;
        for (LGraph component : components) {
            KVector componentSize = component.getSize();
            maxRowWidth = Math.max(maxRowWidth, componentSize.x);
            totalArea += componentSize.x * componentSize.y;
        }
        
        maxRowWidth = Math.max(maxRowWidth, (float) Math.sqrt(totalArea)
                * components.iterator().next().getProperty(LayeredOptions.ASPECT_RATIO));
        
        // Place nodes iteratively into rows
        double xpos = 0, ypos = 0, highestBox = 0, broadestRow = spacing;
        for (LGraph graph : components) {
            KVector size = graph.getSize();
            
            if (xpos + size.x > maxRowWidth) {
                // Place the graph into the next row
                xpos = 0;
                ypos += highestBox + spacing;
                highestBox = 0;
            }
            
            offsetGraph(graph, xpos, ypos);
            
            broadestRow = Math.max(broadestRow, xpos + size.x);
            highestBox = Math.max(highestBox, size.y);
            
            xpos += size.x + spacing;
        }
        
        return new KVector(broadestRow + spacing, ypos + highestBox + spacing);
    }

}
