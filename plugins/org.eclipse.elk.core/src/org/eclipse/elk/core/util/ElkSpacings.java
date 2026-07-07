/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 ********************************************************************************/
package org.eclipse.elk.core.util;

import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.core.LayoutConfigurator.IPropertyHolderOptionFilter;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.math.DoubleMath;

/**
 * Utility class to conveniently configure core spacing values ({@link CoreOptions#SPACING_*}) for a graph element or a
 * whole graph.
 * 
 * Note that individual layout algorithms may specify additional spacing values and may hence provide their own utility
 * class to configure those spacings. For layered this is definitely the case ({@link LayeredSpacings}).
 */
public final class ElkSpacings {

    private ElkSpacings() { }

    /**
     * Start constructing a new spacing configurator with the passed {@code baseSpacing} as base value for all other
     * spacing options. By default the values for other spacing options {@code o} will be determined based on the ratio
     * between {@link ElkCoreSpacingsBuilder#BASE_SPACING_OPTION}'s default value and {@code o}'s default value.
     */
    public static ElkCoreSpacingsBuilder withBaseValue(final double baseSpacing) {
        return new ElkCoreSpacingsBuilder(baseSpacing);
    }
    
    /**
     * Spacing configuration builder implementation for spacing options defined in {@link CoreOptions}.
     */
    public static final class ElkCoreSpacingsBuilder extends AbstractSpacingsBuilder<ElkCoreSpacingsBuilder> {

        /** See {@link AbstractSpacingsBuilder#baseSpacing} for details. */
        public static final IProperty<Double> BASE_SPACING_OPTION = CoreOptions.SPACING_NODE_NODE;
        
        // @formatter:off
        private static final List<IProperty<Double>> DEPENDENT_SPACING_OPTIONS = Lists.newArrayList(
            // Sorted alphabetically
            CoreOptions.SPACING_COMPONENT_COMPONENT,
            CoreOptions.SPACING_EDGE_EDGE,
            CoreOptions.SPACING_EDGE_LABEL,
            CoreOptions.SPACING_EDGE_NODE,
            CoreOptions.SPACING_LABEL_LABEL,
            CoreOptions.SPACING_LABEL_NODE,
            CoreOptions.SPACING_LABEL_PORT_HORIZONTAL,
            CoreOptions.SPACING_LABEL_PORT_VERTICAL,
            CoreOptions.SPACING_NODE_SELF_LOOP,
            CoreOptions.SPACING_PORT_PORT
        );
        // @formatter:on

        static {
            assert BASE_SPACING_OPTION.getDefault() != null : "Base spacing default value must be non-null.";
            // Avoid division by zero.
            assert !DoubleMath.fuzzyEquals(0.0d, BASE_SPACING_OPTION.getDefault(),
                    DOUBLE_EQ_EPSILON) : "Base spacing default value must be different from 0.0d.";
        }
        
        private ElkCoreSpacingsBuilder(final double baseSpacing) {
            super(baseSpacing);
        }

        @Override
        protected ElkCoreSpacingsBuilder thisT() {
            return this;
        }
        
        @Override
        protected IProperty<Double> getBaseSpacingOption() {
            return BASE_SPACING_OPTION;
        }
        
        @Override
        protected Iterable<IProperty<Double>> getDependendSpacingOptions() {
            return DEPENDENT_SPACING_OPTIONS;
        }
    }

    /**
     * An abstract builder to be extended to cater for the needs of individual layout algorithms when it comes to layout
     * options only they are aware of.
     * 
     * @param <T>
     *            The class of the inheriting builder itself.
     */
    public abstract static class AbstractSpacingsBuilder<T extends AbstractSpacingsBuilder<T>> {
     
        /**
         * Delegates to {@link LayoutConfigurator#OPTION_TARGET_FILTER} if the holder is an {@link ElkGraphElement},
         * otherwise returns true.
         * 
         * Required here since the base spacing builder class should not restrict itself to {@link ElkGraphElement}s
         * only and thus allows configuration for arbitrary {@link IPropertyHolder}s.
         */
        public static final IPropertyHolderOptionFilter ELK_OPTION_TARGET_FILTER = (holder, property) -> {
            if (property instanceof ElkGraphElement) {
                return LayoutConfigurator.OPTION_TARGET_FILTER.accept((ElkGraphElement) holder, property);
            } else {
                return true;
            }
        };

        /** Constant for fuzzy double comparisons. */
        protected static final double DOUBLE_EQ_EPSILON = 10e-5;
        

        // Internal builder configuration elements
        
        /** The spacing layout option that shall serve as reference when deriving spacing values for all other spacing
         * layout options. In other words the ratio between this option's default value and a user-specified base
         * value for a "general spacings" yields a factor to be used for the other spacing layout option default values.
         */
        private final double baseSpacing;
        /** Whether the resulting configurator shall overwrite existing spacing option values. */
        private boolean overwrite = false;
        /** Filters to be applied by the resulting configurator. Elements specified directly here are always applied. */
        private List<IPropertyHolderOptionFilter> filters = Lists.newArrayList(ELK_OPTION_TARGET_FILTER);
        
        /**
         * Internal mapping of dependent layout options to factors that must be applied to those option's default 
         * values when setting the spacing value for a concrete node.
         */
        private Map<IProperty<Double>, Double> factorMap = Maps.newHashMap();
                
        /**
         * Constructs a new builder and computes the initial factors based on the passed based value.
         */
        protected AbstractSpacingsBuilder(final double theBaseSpacing) {
            this.baseSpacing = theBaseSpacing;

            // Put the base spacing option in the map as well to allow iterating all options by using this map
            factorMap.put(getBaseSpacingOption(), 1.0d);

            // Initially, and by default, set all factors to the ratio given by default values set in the Core.melk
            getDependendSpacingOptions().forEach(option -> {
                double factor = (getBaseSpacingOption().getDefault() != null && option.getDefault() != null)
                        ? option.getDefault() / getBaseSpacingOption().getDefault()
                        : 1.0d;
                factorMap.put(option, factor);
            });
        }
        
        // Protected

        /**
         * See {@link #baseSpacing}'s javadoc for further details.
         * 
         * @return a spacing option that shall serve as the basis for the configuration, i.e. it will serve as basis to
         *         compute factors to be applied to other spacing values.
         */
        protected abstract IProperty<Double> getBaseSpacingOption();

        /**
         * @return a list of spacing options that shall be affected by the specified base value.
         */
        protected abstract Iterable<IProperty<Double>> getDependendSpacingOptions();
        
        /**
         * @return the filters that shall be applied during configuration.
         */
        protected List<IPropertyHolderOptionFilter> getFilters() {
            return filters;
        }

        /**
         * @return type of the implementing class - necessary to make suppress a warning whenever a method in this base
         *         class returns {@code T}.
         */
        protected abstract T thisT();
        
        /**
         * @return the factor to be used for {@code value} relative to the specified {@link #baseSpacing}.
         */
        protected double computeFactor(final double value) {
            return value / baseSpacing;
        }
        
        // Public API
        
        /**
         * @return an immutable list of the currently configured factors.
         */
        public Map<IProperty<Double>, Double> getFactors() {
            return ImmutableMap.copyOf(factorMap);
        }

        /**
         * Configure the {@code factor} to use for {@code spacingOption} relative to the configured base value.
         * 
         * @return this instance.
         */
        public T withFactor(final IProperty<Double> spacingOption, final double factor) {
            if (!factorMap.keySet().contains(spacingOption)) {
                throw new IllegalArgumentException(
                        "'" + spacingOption.getId() + "' is not a configurable spacing option.");
            }
            if (factor < 0.0) {
                throw new IllegalArgumentException(
                        "The factor for '" + spacingOption.getId() + "' must not be negative (" + factor + ").");
            }
            if (spacingOption.getId().equals(getBaseSpacingOption().getId())) {
                throw new IllegalArgumentException("'" + spacingOption.getId()
                        + "' is the base spacing option not allowed to use with 'withFactor'.");
            }
            factorMap.put(spacingOption, factor);
            return thisT();
        }

        /**
         * Configure the {@code value} to use for {@code spacingOption} relative to the configured base value.
         * 
         * @return this instance.
         */
        public T withValue(final IProperty<Double> spacingOption, final double value) {
            if (!factorMap.keySet().contains(spacingOption)) {
                throw new IllegalArgumentException(
                        "'" + spacingOption.getId() + "' is not a configurable spacing option.");
            }
            if (value < 0.0) {
                throw new IllegalArgumentException(
                        "The value for '" + spacingOption.getId() + "' must not be negative (" + value + ").");
            }
            if (spacingOption.getId().equals(getBaseSpacingOption().getId())) {
                throw new IllegalArgumentException("'" + spacingOption.getId()
                        + "' is the base spacing option not allowed to use with 'withValue'.");
            }
            factorMap.put(spacingOption, computeFactor(value));
            return thisT();
        }
        
        /**
         * If {@code overwrite} is set to false, existing spacing option values will not be overwritten when configuring
         * an element.
         * 
         * @return this instance.
         */
        public T withOverwrite(final boolean shallOverwrite) {
            this.overwrite = shallOverwrite;
            return thisT();
        }
        
        /**
         * @return a function accepting a {@link IPropertyHolder} and applying the assembled spacing configuration to it
         *         when called.
         */
        public Consumer<IPropertyHolder> build() {
            if (!this.overwrite) {
                filters.add(LayoutConfigurator.NO_OVERWRITE_HOLDER);
            }
            
            // Early exit if we are not allowed to overwrite any options and if the specified base matches the default
            // values anyway
            if (!overwrite
                    && DoubleMath.fuzzyEquals(baseSpacing, getBaseSpacingOption().getDefault(), DOUBLE_EQ_EPSILON)) {
                return e -> { };
            }

            return element -> {
                // Set each spacing option's value based on the computed factors. Remember that the base value option 
                // is included in the map. Only apply to graph elements:
                //  - that actually support this option
                //  - for which the option has not been set so far
                factorMap.keySet().stream()
                    .filter(p -> filters.stream().allMatch(filter -> filter.accept(element, p)))
                    .forEach(p -> {
                        element.setProperty(p, factorMap.get(p) * baseSpacing);
                    });
            };
        }
        
        /**
         * Apply the built configuration to the passed {@code holder}. Note that no hierarchy is traversed. Use an
         * {@link IGraphElementVisitor} returned by {@link #toVisitor()} for that.
         * 
         * @param holder
         */
        public void apply(final IPropertyHolder holder) {
           build().accept(holder);
        }
        
        /**
         * @return a new {@link IGraphElementVisitor} that applies the built configuration to each node it visits.
         */
        public IGraphElementVisitor toVisitor() {
            Consumer<IPropertyHolder> consumer = build();
            return e -> consumer.accept(e);
        }
     
    }

}
