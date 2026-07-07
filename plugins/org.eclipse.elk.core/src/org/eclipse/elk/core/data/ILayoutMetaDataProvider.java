/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.data;

/**
 * Provider for meta data of layout algorithms and layout options.
 */
public interface ILayoutMetaDataProvider {
    
    /**
     * Apply this provider by registering all contained meta data into the given registry instance.
     */
    void apply(Registry registry);
    
    /**
     * Registry for layout meta data.
     */
    interface Registry {
        
        /**
         * Register a layout algorithm.
         */
        void register(LayoutAlgorithmData algorithmData);
        
        /**
         * Register a layout option.
         */
        void register(LayoutOptionData optionData);
        
        /**
         * Register a layout category.
         */
        void register(LayoutCategoryData categoryData);
        
        /**
         * Specify a dependency between two layout options. The source option must be made visible in the UI
         * only if the target option is set to the specified value. If no value is given, the requirement for
         * the target option is to be non-null and non-empty.
         * 
         * @param sourceOption the dependent layout option id
         * @param targetOption the option id for which to check the required value
         * @param requiredValue a value to require for the target option, or {@code null}
         */
        void addDependency(String sourceOption, String targetOption, Object requiredValue);
        
        /**
         * Specify support of a layout algorithm for the given layout option.
         * 
         * @param algorithm a layout algorithm id
         * @param option the option id for which to add support
         * @param defaultValue the default value applied by the algorithm, or {@code null}
         */
        void addOptionSupport(String algorithm, String option, Object defaultValue);
        
    }

}
