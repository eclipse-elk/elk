/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
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

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.core.options.PortSide;

/**
 * Model order component group that saves additionally to the port side the order of the components.
 *
 */
public class ModelOrderComponentGroup extends ComponentGroup {
    
    /**
     * Update constraints map with all constraints that occur if the key is inserted in a component group with the 
     * value already in it. This contains only additional sides to the CONSTRAINTS of the super class.
     */
    protected static final Map<Set<PortSide>, Set<Set<PortSide>>> MODEL_ORDER_CONSTRAINTS = new HashMap<>();

    private static void putModelOrderConstraint(final Set<PortSide> key, final Set<PortSide> value) {
        MODEL_ORDER_CONSTRAINTS.computeIfAbsent(key, k -> new HashSet<>()).add(value);
    }

    static {
        // Key is inserted in component group with value
        // E.g. NORTH can not be put in the same component group as NONE since it would be placed above NONE, which
        // is against the model order.
        putModelOrderConstraint(SIDES_NORTH, SIDES_NONE);
        putModelOrderConstraint(SIDES_WEST, SIDES_NONE);
        putModelOrderConstraint(SIDES_NORTH_EAST, SIDES_NONE);
        putModelOrderConstraint(SIDES_NORTH_WEST, SIDES_NONE);
        putModelOrderConstraint(SIDES_NORTH_SOUTH_WEST, SIDES_NONE);
        putModelOrderConstraint(SIDES_NORTH_EAST_WEST, SIDES_NONE);
        
        putModelOrderConstraint(SIDES_NORTH_WEST, SIDES_NORTH);
        
        putModelOrderConstraint(SIDES_NONE, SIDES_EAST);
        putModelOrderConstraint(SIDES_NORTH, SIDES_EAST);
        putModelOrderConstraint(SIDES_WEST, SIDES_EAST);
        putModelOrderConstraint(SIDES_NORTH_EAST, SIDES_EAST);
        putModelOrderConstraint(SIDES_NORTH_SOUTH, SIDES_EAST);
        putModelOrderConstraint(SIDES_NORTH_WEST, SIDES_EAST);
        putModelOrderConstraint(SIDES_NORTH_SOUTH_WEST, SIDES_EAST);
        putModelOrderConstraint(SIDES_NORTH_EAST_WEST, SIDES_EAST);
        putModelOrderConstraint(SIDES_EAST_WEST, SIDES_EAST); // XXX
        
        putModelOrderConstraint(SIDES_NONE, SIDES_SOUTH);
        putModelOrderConstraint(SIDES_NORTH, SIDES_SOUTH);
        putModelOrderConstraint(SIDES_EAST, SIDES_SOUTH);
        putModelOrderConstraint(SIDES_WEST, SIDES_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_EAST, SIDES_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_SOUTH, SIDES_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_WEST, SIDES_SOUTH);
        putModelOrderConstraint(SIDES_EAST_WEST, SIDES_SOUTH);
        putModelOrderConstraint(SIDES_SOUTH_WEST, SIDES_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_SOUTH_WEST, SIDES_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_EAST_SOUTH, SIDES_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_EAST_WEST, SIDES_SOUTH);
        
        putModelOrderConstraint(SIDES_NORTH, SIDES_WEST);
        putModelOrderConstraint(SIDES_NORTH_EAST, SIDES_WEST);
        putModelOrderConstraint(SIDES_NORTH_WEST, SIDES_WEST);
        putModelOrderConstraint(SIDES_NORTH_EAST_WEST, SIDES_WEST);
        
        putModelOrderConstraint(SIDES_NORTH, SIDES_NORTH_EAST);
        putModelOrderConstraint(SIDES_WEST, SIDES_NORTH_EAST);
        putModelOrderConstraint(SIDES_NORTH_WEST, SIDES_NORTH_EAST);
        putModelOrderConstraint(SIDES_NORTH_EAST, SIDES_NORTH_EAST);
        putModelOrderConstraint(SIDES_NORTH_SOUTH_WEST, SIDES_NORTH_EAST);
        
        // NW has nothing since it is in the first slot

        // Only conflicts since it is in the last slot
        putModelOrderConstraint(SIDES_NONE, SIDES_EAST_SOUTH);
        putModelOrderConstraint(SIDES_NORTH, SIDES_EAST_SOUTH);
        putModelOrderConstraint(SIDES_EAST, SIDES_EAST_SOUTH);
        putModelOrderConstraint(SIDES_SOUTH, SIDES_EAST_SOUTH);
        putModelOrderConstraint(SIDES_WEST, SIDES_EAST_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_EAST, SIDES_EAST_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_SOUTH, SIDES_EAST_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_WEST, SIDES_EAST_SOUTH);
        putModelOrderConstraint(SIDES_SOUTH_WEST, SIDES_EAST_SOUTH);
//        putModelOrderConstraint(SIDES_EAST_SOUTH, SIDES_EAST_SOUTH);
        putModelOrderConstraint(SIDES_EAST_WEST, SIDES_EAST_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_EAST_WEST, SIDES_EAST_SOUTH);
//        putModelOrderConstraint(SIDES_NORTH_EAST_SOUTH, SIDES_EAST_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_SOUTH_WEST, SIDES_EAST_SOUTH);
//        putModelOrderConstraint(SIDES_EAST_SOUTH_WEST, SIDES_EAST_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_EAST_SOUTH);

        putModelOrderConstraint(SIDES_NONE, SIDES_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH, SIDES_SOUTH_WEST);
        putModelOrderConstraint(SIDES_EAST, SIDES_SOUTH_WEST);
        putModelOrderConstraint(SIDES_WEST, SIDES_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH_EAST, SIDES_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH_SOUTH, SIDES_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH_WEST, SIDES_SOUTH_WEST);
        putModelOrderConstraint(SIDES_EAST_WEST, SIDES_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH_EAST_WEST, SIDES_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH_EAST_SOUTH, SIDES_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_SOUTH_WEST);
        
        putModelOrderConstraint(SIDES_NORTH, SIDES_EAST_WEST);
        putModelOrderConstraint(SIDES_WEST, SIDES_EAST_WEST);
        putModelOrderConstraint(SIDES_NORTH_EAST, SIDES_EAST_WEST);
        putModelOrderConstraint(SIDES_NORTH_WEST, SIDES_EAST_WEST);
        putModelOrderConstraint(SIDES_SOUTH_WEST, SIDES_EAST_WEST);
        putModelOrderConstraint(SIDES_NORTH_EAST_WEST, SIDES_EAST_WEST);
        putModelOrderConstraint(SIDES_NORTH_SOUTH_WEST, SIDES_EAST_WEST);

        // NEW no additional conflicts
        
        putModelOrderConstraint(SIDES_NONE, SIDES_EAST_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH, SIDES_EAST_SOUTH_WEST);
        putModelOrderConstraint(SIDES_EAST, SIDES_EAST_SOUTH_WEST);
        putModelOrderConstraint(SIDES_WEST, SIDES_EAST_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH_EAST, SIDES_EAST_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH_SOUTH, SIDES_EAST_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH_WEST, SIDES_EAST_SOUTH_WEST);
        putModelOrderConstraint(SIDES_EAST_WEST, SIDES_EAST_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH_EAST_WEST, SIDES_EAST_SOUTH_WEST);
        
        putModelOrderConstraint(SIDES_NORTH, SIDES_NORTH_SOUTH_WEST);
        putModelOrderConstraint(SIDES_EAST, SIDES_NORTH_SOUTH_WEST);
        putModelOrderConstraint(SIDES_SOUTH, SIDES_NORTH_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH_EAST, SIDES_NORTH_SOUTH_WEST);
        
        putModelOrderConstraint(SIDES_NONE, SIDES_NORTH_EAST_SOUTH);
        putModelOrderConstraint(SIDES_NORTH, SIDES_NORTH_EAST_SOUTH);
        putModelOrderConstraint(SIDES_SOUTH, SIDES_NORTH_EAST_SOUTH);
        putModelOrderConstraint(SIDES_WEST, SIDES_NORTH_EAST_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_EAST, SIDES_NORTH_EAST_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_SOUTH, SIDES_NORTH_EAST_SOUTH);
        putModelOrderConstraint(SIDES_NORTH_WEST, SIDES_NORTH_EAST_SOUTH);

        putModelOrderConstraint(SIDES_NORTH_WEST, SIDES_NORTH_EAST_SOUTH_WEST);
        putModelOrderConstraint(SIDES_NORTH_EAST, SIDES_NORTH_EAST_SOUTH_WEST);
        
        // Conflicts that seem solvable but that arise since the order of C, EW, W, E is fix
        // EW is before C and W and E
        putModelOrderConstraint(SIDES_EAST_WEST, SIDES_NONE);
        putModelOrderConstraint(SIDES_EAST_WEST, SIDES_WEST);
        putModelOrderConstraint(SIDES_EAST_WEST, SIDES_EAST);
        
        // Conflicts that seem solvable but that arise since the order of C, NS, N, S is fix
        putModelOrderConstraint(SIDES_NORTH_SOUTH, SIDES_NONE);
        putModelOrderConstraint(SIDES_NORTH_SOUTH, SIDES_NORTH);
        putModelOrderConstraint(SIDES_NORTH_SOUTH, SIDES_SOUTH);
        
    }
    
    
    
    private List<LGraph> componentOrder = new LinkedList<>();
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
            components
                .computeIfAbsent(component.getProperty(InternalProperties.EXT_PORT_CONNECTIONS),
                        k -> new ArrayList<>())
                .add(component);
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
        Set<Set<PortSide>> constraints = CONSTRAINTS.getOrDefault(candidateSides, Collections.emptySet());
        Set<Set<PortSide>> modelOrderConstraints = MODEL_ORDER_CONSTRAINTS.getOrDefault(candidateSides,
                Collections.emptySet());

        for (Set<PortSide> constraint : constraints) {
            if (!components.getOrDefault(constraint, Collections.emptyList()).isEmpty()) {
                // A component with a conflicting external port side combination exists
                return false;
            }
        }
        for (Set<PortSide> constraint : modelOrderConstraints) {
            if (!components.getOrDefault(constraint, Collections.emptyList()).isEmpty()) {
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
