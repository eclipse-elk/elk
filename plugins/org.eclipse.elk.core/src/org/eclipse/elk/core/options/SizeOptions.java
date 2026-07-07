/*******************************************************************************
 * Copyright (c) 2012, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.options;

/**
 * Size constraint options that modify how size constraints of a node are applied. Interpreting a
 * size option set on a node only makes sense when its size constraint is taken into consideration
 * as well. <b>This enumeration is not set directly on {@link CoreOptions#SIZE_OPTIONS}; instead,
 * an {@code EnumSet} over this enumeration is used there.</b>
 * 
 * <p><i>Note:</i> Layout algorithms may only support a subset of these options.</p>
 *
 * @author cds
 */
public enum SizeOptions {
    
    /**
     * If no minimum size is set on an element, the minimum size options are assumed to be some
     * default value determined by the particular layout algorithm. This option only makes sense
     * if the {@link SizeConstraint#MINIMUM_SIZE} constraint is set.
     */
    DEFAULT_MINIMUM_SIZE,
    
    /**
     * If this option is set and paddings are computed by the algorithm, the minimum size plus the
     * computed padding are a lower bound on the node size. If this option is not set, the minimum size
     * will be applied to the node's whole size regardless of any computed padding. Note that,
     * depending on the algorithm, this option may only apply to non-hierarchical nodes. This option
     * only makes sense if the {@link SizeConstraint#MINIMUM_SIZE} constraint is set.
     */
    MINIMUM_SIZE_ACCOUNTS_FOR_PADDING,
    
    /**
     * With this option set, the padding of nodes will be computed and returned as part of the
     * algorithm's result. If port labels or node labels are placed, they may influence the size of
     * the padding. Note that, depending on the algorithm, this option may only apply to
     * non-hierarchical nodes. This option is independent of the size constraint set on a node.
     */
    COMPUTE_PADDING,
    
    /**
     * If node labels influence the node size, but outside node labels are allowed to overhang, only inside node labels
     * actually influence node size.
     */
    OUTSIDE_NODE_LABELS_OVERHANG,
    
    /**
     * By default, ports only use the space available to them, even if that means violating port spacing settings. If
     * this option is active, port spacings are adhered to, even if that means ports extend beyond node boundaries.
     */
    PORTS_OVERHANG,
    
    /**
     * If port labels are taken into consideration, differently sized labels can result in a different amount of space
     * between different pairs of ports. This option causes all ports to be evenly spaced by enlarging the space
     * between every pair of ports to the largest amount of space between any pair of ports.
     */
    UNIFORM_PORT_SPACING,

    /**
     * @deprecated Use {@link PortLabelPlacement#SPACE_EFFICIENT}.
     */
    @Deprecated
    SPACE_EFFICIENT_PORT_LABELS,

    /**
     * By default, inside node labels will be laid out in three rows of three cells, with no relation between the
     * width of cells in different rows. If this option is enabled, the cells will be treated as cells of a table,
     * with equal columns across all rows. This usually results in larger nodes.
     */
    FORCE_TABULAR_NODE_LABELS,
    
    /**
     * If this option is set, the node sizing and label placement code will not make an attempt to achieve a symmetrical
     * layout. With this option inactive, for example, the space reserved for left inside port labels will be the same
     * as for right inside port labels, which would not be the case otherwise. Deactivating this option will also ensure
     * that center node labels will actually be placed in the center.
     */
    ASYMMETRICAL;
    
    /**
     * Returns the enumeration value related to the given ordinal.
     * 
     * @param i ordinal value
     * @return the related enumeration value
     */
    public static SizeOptions valueOf(final int i) {
        return values()[i];
    }
}
