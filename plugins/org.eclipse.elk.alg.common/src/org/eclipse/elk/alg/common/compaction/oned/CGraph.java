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

import java.util.EnumSet;
import java.util.List;

import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.Lists;

/**
 * Representation of a constraint graph that can be compacted using the {@link OneDimensionalCompactor}.
 */
public final class CGraph {
    
    // Variables are public for convenience reasons
    // SUPPRESS CHECKSTYLE NEXT 14 VisibilityModifier
    /** the list of {@link CNode}s modeling the constraints in this graph. */
    public List<CNode> cNodes = Lists.newArrayList();
    /** groups of elements that are supposed to stay in the configuration they are. */
    public List<CGroup> cGroups = Lists.newArrayList();
    /** the directions that are supported for compaction. */
    private EnumSet<Direction> supportedDirections;
    
    /** static constraints between pairs of {@link CNode}s that are to be added by the {@link OneDimensionalCompactor}.
     *  Depending on the compaction direction, the constraints may be applied in reverse. */
    public List<Pair<CNode, CNode>> predefinedHorizontalConstraints = Lists.newArrayList();
    /** static constraints between pairs of {@link CNode}s that are to be added by the {@link OneDimensionalCompactor}.
     *  Depending on the compaction direction, the constraints may be applied in reverse. */
    public List<Pair<CNode, CNode>> predefinedVerticalConstraints = Lists.newArrayList();
    
    /**
     * Constructor sets the supported directions.
     * 
     * @param supportedDirections
     *          the directions that are supported for compaction
     */
    public CGraph(final EnumSet<Direction> supportedDirections) {
        this.supportedDirections = supportedDirections;
    }
    
    /**
     * If the {@link CGraph} supports compaction in the direction specified by the parameter.
     * @param direction
     *          the direction to check
     * @return if compaction is supported
     */
    public boolean supports(final Direction direction) {
        return supportedDirections.contains(direction);
    }
}
