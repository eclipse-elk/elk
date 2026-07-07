/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.compaction.oned;

import org.eclipse.elk.core.options.Direction;

/**
 * Internal class representing a 4-tuple that, in one application as a 'compaction lock', states for
 * a {@link CNode} if the compaction should be locked in a particular direction.
 */
public final class Quadruplet {
    
    // SUPPRESS CHECKSTYLE NEXT 2 VisibilityModifier
    /** Locking values. */
    public boolean left, right, up, down;
    
    /**
     * The lock defaults to false.
     */
    public Quadruplet() {
        set(false, false, false, false);
    }
    
    /**
     * This constructor initializes the lock.
     * 
     * @param l left
     * @param r right
     * @param u up
     * @param d down
     */
    public Quadruplet(final boolean l, final boolean r, final boolean u, final boolean d) {
        set(l, r, u, d);
    }
    
    /**
     * Sets the lock.
     * 
     * @param l left
     * @param r right
     * @param u up
     * @param d down
     */
    public void set(final boolean l, final boolean r, final boolean u, final boolean d) {
        left = l;
        right = r;
        up = u;
        down = d;
    }
    
    /**
     * Modifies this {@link Quadruplet} to represent the element-wise <em>or</em> 
     * with {@code other}.
     * 
     * @param other
     */
    public void applyOr(final Quadruplet other) {
       this.left |= other.left;
       this.right |= other.right;
       this.up |= other.up;
       this.down |= other.down;
    }
    
    /**
     * Sets the lock in a specific {@link Direction}.
     * 
     * @param value 
     *          the desired state
     * @param direction the {@link Direction} of compaction
     */
    public void set(final boolean value, final Direction direction) {
        switch (direction) {
        case LEFT:
            left = value;
            break;
            
        case RIGHT:
            right = value;
            break;
            
        case UP:
            up = value;
            break;
            
        case DOWN:
            down = value;
            break;

        default:
            break;
        }
    }
    
    /**
     * Returns the state for a {@link Direction}.
     * 
     * @param direction
     *          the {@link Direction}
     * @return the state
     */
    public boolean get(final Direction direction) {
        switch (direction) {
        case LEFT:
            return left;
            
        case RIGHT:
            return right;
            
        case UP:
            return up;
            
        case DOWN:
            return down;

        default:
            return false;
        }
    }
}
