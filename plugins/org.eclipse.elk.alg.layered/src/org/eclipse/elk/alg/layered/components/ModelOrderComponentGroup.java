/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
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
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.core.options.PortSide;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

/**
 * Model order component group that saves additionally to the port side the order of the components.
 *
 */
public class ModelOrderComponentGroup extends ComponentGroup {
    
    /**
     * Update constraints map with all constraints that occur if the key is inserted in a component group with the 
     * value already in it. This contains only additional sides to the CONSTRAINTS of the super class.
     */
    protected static final Multimap<Set<PortSide>, Set<PortSide>> MODEL_ORDER_CONSTRAINTS = HashMultimap.create();
    
    static {
        // Key is inserted in component group with value
        // E.g. NORTH can not be put in the same component group as NONE since it would be placed above NONE, which
        // is against the model order.
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH, SIDES_NONE);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_WEST, SIDES_NONE);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_NONE);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_NONE);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_NONE);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_NONE);
        
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_NORTH);
        
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NONE, SIDES_EAST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH, SIDES_EAST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_WEST, SIDES_EAST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_EAST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_EAST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_EAST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_EAST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_EAST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST_WEST, SIDES_EAST); // XXX
        
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NONE, SIDES_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH, SIDES_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST, SIDES_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_WEST, SIDES_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST_WEST, SIDES_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_SOUTH_WEST, SIDES_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH, SIDES_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_SOUTH);
        
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH, SIDES_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_WEST);
        
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH, SIDES_NORTH_EAST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_WEST, SIDES_NORTH_EAST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_NORTH_EAST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_NORTH_EAST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_NORTH_EAST);
        
        // NW has nothing since it is in the first slot

        // Only conflicts since it is in the last slot
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NONE, SIDES_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH, SIDES_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST, SIDES_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_SOUTH, SIDES_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_WEST, SIDES_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_SOUTH_WEST, SIDES_EAST_SOUTH);
//        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST_SOUTH, SIDES_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST_WEST, SIDES_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_EAST_SOUTH);
//        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH, SIDES_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_EAST_SOUTH);
//        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST_SOUTH_WEST, SIDES_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_EAST_SOUTH);

        MODEL_ORDER_CONSTRAINTS.put(SIDES_NONE, SIDES_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH, SIDES_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST, SIDES_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_WEST, SIDES_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST_WEST, SIDES_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH, SIDES_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_SOUTH_WEST);
        
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH, SIDES_EAST_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_WEST, SIDES_EAST_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_EAST_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_EAST_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_SOUTH_WEST, SIDES_EAST_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_EAST_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_EAST_WEST);

        // NEW no additional conflicts
        
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NONE, SIDES_EAST_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH, SIDES_EAST_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST, SIDES_EAST_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_WEST, SIDES_EAST_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_EAST_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_EAST_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_EAST_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST_WEST, SIDES_EAST_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_EAST_SOUTH_WEST);
        
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH, SIDES_NORTH_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST, SIDES_NORTH_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_SOUTH, SIDES_NORTH_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_NORTH_SOUTH_WEST);
        
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NONE, SIDES_NORTH_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH, SIDES_NORTH_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_SOUTH, SIDES_NORTH_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_WEST, SIDES_NORTH_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_NORTH_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_NORTH_EAST_SOUTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_NORTH_EAST_SOUTH);

        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_NORTH_EAST_SOUTH_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_NORTH_EAST_SOUTH_WEST);
        
        // Conflicts that seem solvable but that arise since the order of C, EW, W, E is fix
        // EW is before C and W and E
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST_WEST, SIDES_NONE);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST_WEST, SIDES_WEST);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_EAST_WEST, SIDES_EAST);
        
        // Conflicts that seem solvable but that arise since the order of C, NS, N, S is fix
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_NONE);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_NORTH);
        MODEL_ORDER_CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_SOUTH);
        
    }
    
    
    
    private List<LGraph> componentOrder = Lists.newLinkedList();
    ///////////////////////////////////////////////////////////////////////////////
    // Constructors
    
    /**
     * Constructs a new, empty component group.
     */
    public ModelOrderComponentGroup() {
        
    }
    
    /**
     * Constructs a new component group with the given initial component. This is equivalent to
     * constructing an empty component group and then calling {@link #add(LGraph)}.
     * 
     * @param component the component to be added to the group.
     */
    public ModelOrderComponentGroup(final LGraph component) {
        add(component);
        getComponentOrder().add(component);
    }

    
    ///////////////////////////////////////////////////////////////////////////////
    // Component Management
    
    /**
     * Tries to add the given component to the group. Before adding the component, a call to
     * {@link #canAdd(LGraph)} determines if the component can actually be added to this
     * group.
     * 
     * @param component the component to be added to this group.
     * @return {@code true} if the component was successfully added, {@code false} otherwise.
     */
    public boolean add(final LGraph component) {
        if (canAdd(component)) {
            components.put(
                component.getProperty(InternalProperties.EXT_PORT_CONNECTIONS),
                component);
            getComponentOrder().add(component);
            return true;
        } else {
            return false;
        }
    }
    
    /**
     * Checks whether this group has enough space left to add a given component.
     * 
     * @param component the component to be added to the group.
     * @return {@code true} if the group has enough space left to add the component, {@code false}
     *         otherwise.
     */
    protected boolean canAdd(final LGraph component) {
        // Check if we have a component with incompatible external port sides
        Set<PortSide> candidateSides = component.getProperty(InternalProperties.EXT_PORT_CONNECTIONS);
        Collection<Set<PortSide>> constraints = CONSTRAINTS.get(candidateSides);
        Collection<Set<PortSide>> modelOrderConstraints = MODEL_ORDER_CONSTRAINTS.get(candidateSides);
        
        for (Set<PortSide> constraint : constraints) {
            if (!components.get(constraint).isEmpty()) {
                // A component with a conflicting external port side combination exists
                return false;
            }
        }
        for (Set<PortSide> constraint : modelOrderConstraints) {
            if (!components.get(constraint).isEmpty()) {
                // A component with a conflicting external port side combination exists considering model order
                return false;
            }
        }
        
        // We haven't found any conflicting components
        return true;
    }

    /**
     * @return the componentOrder
     */
    public List<LGraph> getComponentOrder() {
        return componentOrder;
    }
}
