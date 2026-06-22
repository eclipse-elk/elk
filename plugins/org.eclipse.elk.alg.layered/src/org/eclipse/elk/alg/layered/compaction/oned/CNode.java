/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.compaction.oned;

import java.util.List;
import java.util.ArrayList;

import org.eclipse.elk.alg.layered.graph.LGraphElement;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;



/**
 * Internal representation of a node in the constraint graph.
 * 
 * For instance, this class is extended to handle specific
 * {@link org.eclipse.elk.alg.layered.graph.LGraphElement LGraphElement}s.
 * 
 * @see org.eclipse.elk.alg.layered.intermediate.compaction.CLNode CLNode
 * @see org.eclipse.elk.alg.layered.intermediate.compaction.CLEdge CLEdge
 */
public abstract class CNode {
    
    // Variables are public for convenience reasons since this class is used internally only.
    // SUPPRESS CHECKSTYLE NEXT 28 VisibilityModifier
    /** containing {@link CGroup}. */
    public CGroup cGroup;
    /** refers to the parent node of a north/south segment. */
    public CNode parentNode = null;
    /** representation of constraints. */
    public List<CNode> constraints = new ArrayList<>();
    /** the area occupied by this element including margins for ports and labels. */
    public ElkRectangle hitbox;
    /** offset to the root position of the containing {@link CGroup} . */
    public KVector cGroupOffset = new KVector();
    /** leftmost possible position for this {@link CNode} to be drawn. 
     *  This position can be intermediate and is increased to its final value by updateStartPos(). */
    public double startPos = Double.NEGATIVE_INFINITY;
    /** flags a {@link CNode} to be repositioned in the case of left/right balanced compaction. */
    public boolean reposition = true;
    /** a 4-tuple stating if the {@link CNode} should locked in a particular direction based on
     *  conditions defined in an extended class. */
    public Quadruplet lock = new Quadruplet();
    /** Whether no spacing should be applied to a certain side of this node. */
    public Quadruplet spacingIgnore = new Quadruplet();
    /** An id for public use. There is no warranty, use at your own risk. */
    public int id;

    /**
     * Returns the required horizontal spacing to the specified {@link CNode}.
     * 
     * @return the spacing
     */
    public abstract double getHorizontalSpacing();
    
    /**
     * Returns the required vertical spacing to the specified {@link CNode}.
     * 
     * @return the spacing
     */
    public abstract double getVerticalSpacing();

    /**
     * Getter for the position.
     * 
     * @return position of the hitbox
     */
    public double getPosition() {
        return hitbox.x;
    }

    /**
     * Applies the compacted starting position to the hitbox. Used after compaction to allow
     * reverse transformation of hitboxes.
     */
    public void applyPosition() {
        hitbox.x = startPos;
    }

    /**
     * Sets the position of the {@link LGraphElement} according to the hitbox.
     */
    public abstract void applyElementPosition();

    /**
     * Returns the position of the {@link LGraphElement}.
     * 
     * @return the position
     */
    public abstract double getElementPosition();
    
    /**
     * @return an svg representation of this {@link CNode} to the output for debugging.
     */
    protected String getDebugSVG() {
        StringBuffer sb = new StringBuffer();
        sb.append("<rect x=\"" + this.hitbox.x + "\" y=\"" + this.hitbox.y + "\" width=\""
                + Math.max(1, this.hitbox.width) + "\" height=\"" + Math.max(1, this.hitbox.height)
                + "\" fill=\"" + (this.reposition ? "green" : "orange")
                + "\" stroke=\"black\" opacity=\"0.5\"/>");
        sb.append("<text x=\"" + (this.hitbox.x + 2) + "\" y=\""
                + (this.hitbox.y + 2 * 2 * 2 + 2 + 1) + "\">" + "(" + Math.round(this.hitbox.x)
                + ", " + Math.round(this.hitbox.y) + ")\n" + this + "</text>");
        return sb.toString();
    }
}
