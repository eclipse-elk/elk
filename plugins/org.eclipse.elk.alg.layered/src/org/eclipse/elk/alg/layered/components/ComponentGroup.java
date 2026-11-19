/*******************************************************************************
 * Copyright (c) 2012, 2015 Kiel University and others.
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
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.core.options.PortSide;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

/**
 * Represents a group of connected components grouped for layout purposes.
 * 
 * <p>A component group is conceptually divided into nine sectors as such: (the nine sectors are
 * enumerated in the {@link ComponentGroupSector} enumeration)</p>
 * <pre>
 *   +----+----+----+
 *   | nw | n  | ne |
 *   +----+----+----+
 *   | w  | c  | e  |
 *   +----+----+----+
 *   | sw | s  | se |
 *   +----+----+----+
 * </pre>
 * <p>The port sides of external ports a component connects to determines which sector(s) it will
 * occupy. This is best illustrated by some examples:</p>
 * <ul>
 *   <li>Let {@code c} be a component connected to a northern port. Then {@code c} would be placed in
 *       the {@code n} sector.</li>
 *   <li>Let {@code c} be a comopnent connected to a southern port and to an eastern port. Then
 *       {@code c} would be placed in the {@code se} sector.</li>
 *   <li>Let {@code c} be a component connected to no port at all. Then {@code c} would be placed in the
 *       {@code c} sector.</li>
 *   <li>Let {@code c} be a component connected to a western and to an eastern port. Then {@code c}
 *       would be placed in the {@code w}, {@code c}, and {@code e} sectors. If {@code c} would also
 *       connected to a southern port, it would also occupy the {@code sw}, {@code sc}, and {@code se}
 *       sectors.</li>
 * </ul>
 * <p>With this placement comes a bunch of constraints. For example, for a component to occupy the
 * top three sectors, none of them must be occupied by another component yet. If the addition of a
 * component to this group would cause a constraint to be violated, it cannot be added.</p>
 * 
 * <p>This class is not supposed to be public, but needs to be for JUnit tests to find it.</p>
 */
public class ComponentGroup {
    
    ///////////////////////////////////////////////////////////////////////////////
    // Constants
    
    // External Port Connection Constraints
    
    /**
     * A map of constraints used to decide whether a component can be placed in this group.
     * 
     * <p>For a new component that is to be placed in this group, the set of external port sides
     * it connects to implies which sets of port sides of other components it is compatible to.
     * For instance, a component connecting to a northern and an eastern port requires that no
     * other component connects to this particular combination of ports. This map maps sets of
     * port sides to a list of port side sets that must not already exist in this group for a
     * component to be added.</p>
     */
    protected static final Multimap<Set<PortSide>, Set<PortSide>> CONSTRAINTS = HashMultimap.create();
    
    static {
        // Setup constraints
        CONSTRAINTS.put(SIDES_NONE, SIDES_NORTH_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_WEST, SIDES_NORTH_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_WEST, SIDES_NORTH_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_EAST, SIDES_NORTH_EAST_SOUTH);
        CONSTRAINTS.put(SIDES_EAST, SIDES_NORTH_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH, SIDES_NORTH_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH, SIDES_NORTH_EAST_WEST);
        CONSTRAINTS.put(SIDES_SOUTH, SIDES_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_SOUTH, SIDES_NORTH_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_EAST_WEST);
        CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_NORTH_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_NORTH_EAST_WEST);
        CONSTRAINTS.put(SIDES_NORTH_SOUTH, SIDES_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_EAST_WEST, SIDES_NORTH_SOUTH);
        CONSTRAINTS.put(SIDES_EAST_WEST, SIDES_NORTH_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_EAST_WEST, SIDES_NORTH_EAST_SOUTH);
        CONSTRAINTS.put(SIDES_EAST_WEST, SIDES_NORTH_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_NORTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_NORTH_EAST_WEST);
        CONSTRAINTS.put(SIDES_NORTH_WEST, SIDES_NORTH_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_NORTH_EAST);
        CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_NORTH_EAST_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST, SIDES_NORTH_EAST_SOUTH);
        CONSTRAINTS.put(SIDES_SOUTH_WEST, SIDES_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_SOUTH_WEST, SIDES_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_SOUTH_WEST, SIDES_NORTH_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_EAST_SOUTH, SIDES_EAST_SOUTH);
        CONSTRAINTS.put(SIDES_EAST_SOUTH, SIDES_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_EAST_SOUTH, SIDES_NORTH_EAST_SOUTH);
        CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_NORTH);
        CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_NORTH_SOUTH);
        CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_NORTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_NORTH_EAST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_NORTH_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_NORTH_EAST_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_NORTH_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_WEST, SIDES_NORTH_EAST_SOUTH);
        CONSTRAINTS.put(SIDES_EAST_SOUTH_WEST, SIDES_SOUTH);
        CONSTRAINTS.put(SIDES_EAST_SOUTH_WEST, SIDES_NORTH_SOUTH);
        CONSTRAINTS.put(SIDES_EAST_SOUTH_WEST, SIDES_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_EAST_SOUTH_WEST, SIDES_EAST_SOUTH);
        CONSTRAINTS.put(SIDES_EAST_SOUTH_WEST, SIDES_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_EAST_SOUTH_WEST, SIDES_NORTH_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_EAST_SOUTH_WEST, SIDES_NORTH_EAST_SOUTH);
        CONSTRAINTS.put(SIDES_EAST_SOUTH_WEST, SIDES_NORTH_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_WEST);
        CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_EAST_WEST);
        CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_NORTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_NORTH_EAST_WEST);
        CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_NORTH_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_SOUTH_WEST, SIDES_NORTH_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH, SIDES_EAST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH, SIDES_EAST_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH, SIDES_NORTH_EAST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH, SIDES_EAST_SOUTH);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH, SIDES_NORTH_EAST_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH, SIDES_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH, SIDES_NORTH_EAST_SOUTH);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH, SIDES_NORTH_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_NONE);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_EAST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_NORTH);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_SOUTH);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_NORTH_SOUTH);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_EAST_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_NORTH_EAST_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_EAST_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_NORTH_SOUTH_WEST);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_NORTH_EAST_SOUTH);
        CONSTRAINTS.put(SIDES_NORTH_EAST_SOUTH_WEST, SIDES_NORTH_EAST_SOUTH_WEST);
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // Variables
    
    /**
     * A map mapping external port side combinations to components in this group.
     */
    protected Multimap<Set<PortSide>, LGraph> components = ArrayListMultimap.create();
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // Constructors
    
    /**
     * Constructs a new, empty component group.
     */
    public ComponentGroup() {
        
    }
    
    /**
     * Constructs a new component group with the given initial component. This is equivalent to
     * constructing an empty component group and then calling {@link #add(LGraph)}.
     * 
     * @param component the component to be added to the group.
     */
    public ComponentGroup(final LGraph component) {
        add(component);
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
        
        for (Set<PortSide> constraint : constraints) {
            if (!components.get(constraint).isEmpty()) {
                // A component with a conflicting external port side combination exists
                return false;
            }
        }
        
        // We haven't found any conflicting components
        return true;
    }
    
    /**
     * Returns all port sides in this component group.
     * 
     * @return all port sides in this component group.
     */
    public Collection<Set<PortSide>> getPortSides() {
        return components.keys();
    }
    
    /**
     * Returns all components in this component group.
     * 
     * @return the components in this component group.
     */
    public Collection<LGraph> getComponents() {
        return components.values();
    }
    
    /**
     * Returns the components in this component group connected to external ports on the given set
     * of port sides.
     * 
     * @param connections external port sides the returned components are to be connected to.
     * @return the collection of components. If there are no components, an empty collection is
     *         returned.
     */
    public Collection<LGraph> getComponents(final Set<PortSide> connections) {
        return components.get(connections);
    }
}
