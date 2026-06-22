/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0  
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

import java.util.List;
import java.util.Set;
import java.util.function.BiFunction;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.ElkSpacings.AbstractSpacingsBuilder;
import org.eclipse.elk.core.util.ElkSpacings.ElkCoreSpacingsBuilder;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;

import org.eclipse.elk.alg.layered.compaction.oned.CompareFuzzy;

/**
 * Utility class to conveniently configure core spacing values ({@link LayeredOptions#SPACING_*}) for a graph element or
 * a whole graph.
 */
public final class LayeredSpacings {
    
    private LayeredSpacings() { }

    /**
     * Start constructing a new spacing configurator with the passed {@code baseSpacing} as base value for all other
     * spacing options. By default the values for other spacing options {@code o} will be determined based on the ratio
     * between {@link ElkCoreSpacingsBuilder#BASE_SPACING_OPTION}'s default value and {@code o}'s default value.
     */
    public static LayeredSpacingsBuilder withBaseValue(final double baseSpacing) {
        return new LayeredSpacingsBuilder(baseSpacing);
    }

    /**
     * Spacing configuration builder implementation for spacing options defined in {@link LayeredOptions}.
     */
    public static final class LayeredSpacingsBuilder extends AbstractSpacingsBuilder<LayeredSpacingsBuilder> {

        /** See {@link AbstractSpacingsBuilder#baseSpacing} for details. */
        public static final IProperty<Double> BASE_SPACING_OPTION = LayeredOptions.SPACING_NODE_NODE;
        
        // @formatter:off
        private static final List<IProperty<Double>> DEPENDENT_SPACING_OPTIONS = java.util.Arrays.asList(
            // Sorted alphabetically
            // First the ones inherited from CoreOptions (with possibly overridden default value)
            LayeredOptions.SPACING_COMPONENT_COMPONENT,
            LayeredOptions.SPACING_EDGE_EDGE,
            LayeredOptions.SPACING_EDGE_LABEL,
            LayeredOptions.SPACING_EDGE_NODE,
            LayeredOptions.SPACING_LABEL_LABEL,
            LayeredOptions.SPACING_LABEL_NODE,
            LayeredOptions.SPACING_LABEL_PORT_HORIZONTAL,
            LayeredOptions.SPACING_LABEL_PORT_VERTICAL,
            LayeredOptions.SPACING_NODE_SELF_LOOP,
            LayeredOptions.SPACING_PORT_PORT,
            // Introduced by LayeredOptions 
            LayeredOptions.SPACING_EDGE_EDGE_BETWEEN_LAYERS,
            LayeredOptions.SPACING_EDGE_NODE_BETWEEN_LAYERS,
            LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS
        );
        // @formatter:on

        static {
            assert BASE_SPACING_OPTION.getDefault() != null : "Base spacing default value must be non-null.";
            // Avoid division by zero.
            assert !CompareFuzzy.fuzzyEquals(0.0d, BASE_SPACING_OPTION.getDefault(),
                    DOUBLE_EQ_EPSILON) : "Base spacing default value must be different from 0.0d.";
        }
        
        private LayeredSpacingsBuilder(final double d) {
            super(d);
        }
        
        @Override
        public LayeredSpacingsBuilder thisT() {
            return this;
        }

        @Override
        protected IProperty<Double> getBaseSpacingOption() {
            return BASE_SPACING_OPTION;
        }
        
        /**
         * Note that the super implementation is not simply amended here but completely replaced to properly use the
         * default values from {@link LayeredOptions} for layout options inherited from {@link CoreOptions}.
         * 
         * {@inheritDoc}
         */
        @Override
        protected List<IProperty<Double>> getDependendSpacingOptions() {
            return DEPENDENT_SPACING_OPTIONS;
        }
    }
    
    /**
     * Mimics the behavior of {@link LayoutConfigurator#OPTION_TARGET_FILTER}, just for the {@link LGraph}.
     */
    public static final BiFunction<IPropertyHolder, IProperty<?>, Boolean> OPTION_TARGET_FILTER = (e, property) -> {
        LayoutOptionData optionData = LayoutMetaDataService.getInstance().getOptionData(property.getId());
        if (optionData != null) {
            Set<LayoutOptionData.Target> targets = optionData.getTargets();
            if (e instanceof LGraph) {
                return targets.contains(LayoutOptionData.Target.NODES)
                        || targets.contains(LayoutOptionData.Target.PARENTS);
            } else if (e instanceof LNode) {
                return targets.contains(LayoutOptionData.Target.NODES);
            } else if (e instanceof LEdge) {
                return targets.contains(LayoutOptionData.Target.EDGES);
            } else if (e instanceof LPort) {
                return targets.contains(LayoutOptionData.Target.PORTS);
            } else if (e instanceof LLabel) {
                return targets.contains(LayoutOptionData.Target.LABELS);
            }
        }
        return true;
    };
}
