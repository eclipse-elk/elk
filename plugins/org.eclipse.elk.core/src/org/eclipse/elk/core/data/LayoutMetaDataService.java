/*******************************************************************************
 * Copyright (c) 2008, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.data;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.IndividualSpacings;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.graph.util.ElkReflect;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

/**
 * Singleton class for access to the ELK layout meta data. This class is used globally to retrieve meta data for
 * automatic layout through ELK.
 */
public final class LayoutMetaDataService {

    /** the layout data service instance, which is created lazily. */
    private static LayoutMetaDataService instance;

    private LayoutMetaDataService() {
        initElkReflect();
    }
    
    /**
     * Uses this classes default class loader to load the singleton instance of the layout data service.
     *
     * @return  the singleton instance
     */
    public static synchronized LayoutMetaDataService getInstance() {
        Object loader = null;
        // elkjs-exclude-start
        loader = LayoutMetaDataService.class.getClassLoader();
        // elkjs-exclude-end
        return getInstance(loader);
    }
    
    /**
     * Returns the singleton instance of the layout data service.
     *
     * @param loader
     *      The class loader object. This is not explicitly a class loader since elkjs cannot handle it.
     *
     * @return the singleton instance
     */
    public static synchronized LayoutMetaDataService getInstance(final Object loader) {
        if (instance == null) {
            instance = new LayoutMetaDataService();
            
            // Be sure to load CoreOptions first
            instance.registerLayoutMetaDataProviders(new CoreOptions());

            // Non of the following is available in GWT
            // elkjs-exclude-start
            
            // Invoke service loading to register all meta data providers automatically (this does not work if we're
            // running on Equinox since this will only find services in the realm of this class's class loader)
            for (ILayoutMetaDataProvider provider : java.util.ServiceLoader.load(ILayoutMetaDataProvider.class,
                    (ClassLoader) loader)) {
                instance.registerLayoutMetaDataProviders(provider);
            }

            // Try to make the ELK service plug-in load our services
            try {
                Class.forName("org.eclipse.elk.core.service.ElkServicePlugin");
            } catch (Exception e) {
                // If the service plug-in is not available, that's no problem; we'll simply use our default factory
            }
            // elkjs-exclude-end
        }

        return instance;
    }

    /**
     * Unload any created instance in order to feed the garbage collector.
     */
    public static synchronized void unload() {
        if (instance != null) {
            for (LayoutAlgorithmData algoData : instance.getAlgorithmData()) {
                algoData.getInstancePool().clear();
            }
            instance = null;
        }
    }
    
    /**
     * Registers all basic data types with {@link ElkReflect}. This method should only be called from the outside if
     * layout algorithms are supposed to be called directly, which we don't recommend.
     */
    @SuppressWarnings({ "rawtypes", "unchecked" })
    public static void initElkReflect() {
        // every class that implements IDataObject 
        // or is used for internal properties must be registered here
        
        // IDataObject
        ElkReflect.register(KVector.class, 
                () -> new KVector(), 
                (v) -> ((KVector) v).clone());
        ElkReflect.register(KVectorChain.class, 
                () -> new KVectorChain(), 
                (vc) -> new KVectorChain((KVectorChain) vc));
        ElkReflect.register(ElkMargin.class, 
                () -> new ElkMargin(),
                (m) -> new ElkMargin((ElkMargin) m));
        ElkReflect.register(ElkPadding.class, 
                () -> new ElkPadding(),
                (p) -> new ElkPadding((ElkPadding) p));
        ElkReflect.register(IndividualSpacings.class, 
                () -> new IndividualSpacings(), 
                (is) -> new IndividualSpacings((IndividualSpacings) is));
        
        // Commonly used classes for internal properties
        ElkReflect.register(ArrayList.class, 
                () -> new ArrayList(),
                (al) -> ((ArrayList) al).clone());
        ElkReflect.register(LinkedList.class, 
                () -> new LinkedList(),
                (ll) -> Lists.newLinkedList((LinkedList) ll));
        ElkReflect.register(HashSet.class, 
                () -> new HashSet(),
                (hs) -> Sets.newHashSet((HashSet) hs));
        ElkReflect.register(LinkedHashSet.class, 
                () -> new LinkedHashSet(),
                (hs) -> Sets.newLinkedHashSet((HashSet) hs));
        ElkReflect.register(TreeSet.class, 
                () -> new TreeSet(),
                (ts) -> Sets.newTreeSet((TreeSet) ts));
    }

    /**
     * Mapping of layout provider identifiers to their data instances.
     */
    private final Map<String, LayoutAlgorithmData> layoutAlgorithmMap = Maps.newLinkedHashMap();
    /**
     * Mapping of layout option identifiers to their data instances.
     */
    private final Map<String, LayoutOptionData> layoutOptionMap = Maps.newLinkedHashMap();
    /**
     * Mapping of legacy layout option identifiers to their data instances. Note that the actual layout option data
     * contain the new identifiers.
     */
    private final Map<String, LayoutOptionData> legacyLayoutOptionMap = Maps.newLinkedHashMap();
    /**
     * Mapping of layout category identifiers to their data instances.
     */
    private final Map<String, LayoutCategoryData> layoutCategoryMap = Maps.newLinkedHashMap();
    /**
     * Additional map of layout algorithm suffixes to data instances.
     */
    private final Map<String, LayoutAlgorithmData> algorithmSuffixMap = Maps.newHashMap();
    /**
     * Additional map of layout option suffixes to data instances. For layout options this include the layout option's
     * group.
     */
    private final Map<String, LayoutOptionData> optionSuffixMap = Maps.newHashMap();

    /**
     * Registers the data provided by the given meta data providers with the meta data service.
     *
     * @param providers
     *            the providers to register.
     */
    public void registerLayoutMetaDataProviders(final ILayoutMetaDataProvider... providers) {
        for (ILayoutMetaDataProvider provider : providers) {
            Registry registry = new Registry();
            provider.apply(registry);
            registry.applyDependencies();
        }
        // As soon as new layout options have been registered, we must clear the suffix cache to ensure unique 
        // suffixes from now onwards
        optionSuffixMap.clear();
    }

    /**
     * Returns the layout algorithm data associated with the given identifier.
     *
     * @param id
     *            layout algorithm identifier
     * @return the corresponding layout algorithm data, or {@code null} if there is no algorithm with the given
     *         identifier
     */
    public LayoutAlgorithmData getAlgorithmData(final String id) {
        return layoutAlgorithmMap.get(id);
    }

    /**
     * Returns a data collection for all registered layout algorithms. The collection is unmodifiable.
     *
     * @return collection of registered layout algorithms
     */
    public Collection<LayoutAlgorithmData> getAlgorithmData() {
        return Collections.unmodifiableCollection(layoutAlgorithmMap.values());
    }

    /**
     * Returns a layout algorithm data that has the given suffix in its identifier.
     *
     * @param suffix
     *            a layout algorithm identifier suffix
     * @return the first layout algorithm data that has the given suffix, or {@code null} if no algorithm has that
     *         suffix or the suffix is not unique (multiple algorithms have it)
     */
    public LayoutAlgorithmData getAlgorithmDataBySuffix(final String suffix) {
        if (suffix == null || suffix.isEmpty()) {
            return null;
        }
        LayoutAlgorithmData data = algorithmSuffixMap.get(suffix);

        if (data == null) {
            // Nothing found? Try to find a suitable suffix
            for (LayoutAlgorithmData d : layoutAlgorithmMap.values()) {
                String id = d.getId();
                if (id.endsWith(suffix) && (suffix.length() == id.length()
                        || id.charAt(id.length() - suffix.length() - 1) == '.')) {
                    if (data != null) {
                        // The suffix is not unique, so we don't support it
                        return null;
                    }
                    data = d;
                }
            }

            if (data != null) {
                algorithmSuffixMap.put(suffix, data);
            }
        }

        return data;
    }

    /**
     * Returns the layout algorithm data associated with the given suffix. If none is found, the algorithm with the
     * given default suffix is returned. If that is not found either, returns {@code null}.
     *
     * @param algorithmId
     *            layout algorithm identifier.
     * @param defaultId
     *            fallback layout algorithm identifier. May be {@code null}.
     * @return the layout algorithm data or {@code null} if neither the requested nor the default layout algorithm
     *         are present.
     */
    public LayoutAlgorithmData getAlgorithmDataBySuffixOrDefault(final String algorithmId,
            final String defaultId) {

        // First, look for the requested algorithm
        if (algorithmId != null && !algorithmId.trim().isEmpty()) {
            LayoutAlgorithmData algorithmIdData = getAlgorithmDataBySuffix(algorithmId);

            if (algorithmIdData != null) {
                return algorithmIdData;
            }
        }

        // We haven't found the requested algorithm, so look for the requested default algorithm, if any
        if (defaultId != null && !defaultId.trim().isEmpty()) {
            LayoutAlgorithmData defaultIdData = getAlgorithmDataBySuffix(defaultId);

            if (defaultIdData != null) {
                return defaultIdData;
            }
        }

        // Right, at this point we simply give up
        return null;
    }

    /**
     * Returns the layout option data associated with the given identifier. The identifier may be a legacy identifier.
     *
     * @param id
     *            layout option identifier
     * @return the corresponding layout option data, or {@code null} if there is no option with the given identifier
     */
    public LayoutOptionData getOptionData(final String id) {
        LayoutOptionData data = layoutOptionMap.get(id);
        return data != null ? data : legacyLayoutOptionMap.get(id);
    }

    /**
     * Returns a data collection for all registered layout options. The collection is unmodifiable.
     *
     * @return collection of registered layout options
     */
    public Collection<LayoutOptionData> getOptionData() {
        return Collections.unmodifiableCollection(layoutOptionMap.values());
    }

    /**
     * Returns a layout option data that has the given suffix in its identifier. The identifier may be a legacy
     * identifier.
     *
     * @param suffix
     *            a layout option identifier suffix
     * @return the first layout option data that has the given suffix, or {@code null} if no option has that suffix
     *         or the suffix is not unique (multiple options have it)
     */
    public LayoutOptionData getOptionDataBySuffix(final String suffix) {
        if (suffix == null || suffix.isEmpty()) {
            return null;
        }
        LayoutOptionData data = optionSuffixMap.get(suffix);

        if (data == null) {
            // Nothing found? Try to find a suitable suffix
            for (LayoutOptionData d : layoutOptionMap.values()) {
                String id = d.getId();
                if (id.endsWith(suffix) && (suffix.length() == id.length()
                        || id.charAt(id.length() - suffix.length() - 1) == '.')) {
                    if (data != null) {
                        // The suffix is not unique, so we don't support it
                        return null;
                    }
                    data = d;
                }
            }

            if (data == null) {
                // Still not lucky? Maybe it's a suffix of a legacy id (we should stop supporting these one day...)
                for (LayoutOptionData d : layoutOptionMap.values()) {
                    String[] legacyIds = d.getLegacyIds();
                    if (legacyIds != null) {
                        for (String id : legacyIds) {
                            if (id.endsWith(suffix) && (suffix.length() == id.length()
                                    || id.charAt(id.length() - suffix.length() - 1) == '.')) {
                                if (data != null) {
                                    // The suffix is not unique, so we don't support it
                                    return null;
                                }
                                data = d;
                            }
                        }
                    }
                }
            }

            if (data != null) {
                optionSuffixMap.put(suffix, data);
            }
        }

        return data;
    }

    /**
     * Returns a list of layout options that are suitable for the given layout algorithm and layout option target. The
     * layout algorithm must know the layout options and the target must be active for each option.
     *
     * @param algorithmData
     *            layout algorithm data
     * @param targetType
     *            type of layout option target
     * @return list of suitable layout options
     */
    public List<LayoutOptionData> getOptionData(final LayoutAlgorithmData algorithmData,
            final LayoutOptionData.Target targetType) {

        List<LayoutOptionData> optionDataList = new LinkedList<LayoutOptionData>();
        for (LayoutOptionData optionData : layoutOptionMap.values()) {
            if (algorithmData.knowsOption(optionData) || CoreOptions.ALGORITHM.equals(optionData)) {
                if (optionData.getTargets().contains(targetType)) {
                    optionDataList.add(optionData);
                }
            }
        }
        return optionDataList;
    }

    /**
     * Returns the data instance of the layout category with given identifier.
     *
     * @param id
     *            identifier of the category
     * @return layout category data instance with given identifier, or {@code null} if the layout category is not
     *         registered
     */
    public LayoutCategoryData getCategoryData(final String id) {
        return layoutCategoryMap.get(id);
    }

    /**
     * Returns a list of layout category identifiers and names. The first string in each entry is the identifier, and
     * the second string is the name.
     *
     * @return a list of all layout categories
     */
    public Collection<LayoutCategoryData> getCategoryData() {
        return Collections.unmodifiableCollection(layoutCategoryMap.values());
    }

    /**
     * Class used to register layout services.
     */
    private class Registry implements ILayoutMetaDataProvider.Registry {

        /**
         * Data class for storing dependencies before they are actually registered.
         */
        private class Triple {
            private String firstId;
            private String secondId;
            private Object value;
        }

        private final List<Triple> optionDependencies = new LinkedList<Triple>();
        private final List<Triple> optionSupport = new LinkedList<Triple>();

        @Override
        public void register(final LayoutAlgorithmData algorithmData) {
            layoutAlgorithmMap.put(algorithmData.getId(), algorithmData);
        }

        @Override
        public void register(final LayoutOptionData optionData) {
            // #1 register fully qualified id
            String id = optionData.getId();
            layoutOptionMap.put(id, optionData);

            // #2 register legacy options
            if (optionData.getLegacyIds() != null) {
                for (String legacyId : optionData.getLegacyIds()) {
                    legacyLayoutOptionMap.put(legacyId, optionData);
                }
            }
        }

        @Override
        public void register(final LayoutCategoryData categoryData) {
            layoutCategoryMap.put(categoryData.getId(), categoryData);
        }

        @Override
        public void addDependency(final String sourceOption, final String targetOption, final Object requiredValue) {
            Triple dep = new Triple();
            dep.firstId = sourceOption;
            dep.secondId = targetOption;
            dep.value = requiredValue;
            optionDependencies.add(dep);
        }

        @Override
        public void addOptionSupport(final String algorithm, final String option, final Object defaultValue) {
            Triple sup = new Triple();
            sup.firstId = algorithm;
            sup.secondId = option;
            sup.value = defaultValue;
            optionSupport.add(sup);
        }

        /**
         * Apply all dependencies that have been delivered by layout meta data providers.
         */
        private void applyDependencies() {
            // Go through all registered algorithms and make sure they are
            for (LayoutAlgorithmData algorithm : layoutAlgorithmMap.values()) {
                String categoryId = algorithm.getCategoryId();
                if (categoryId == null) {
                    categoryId = "";
                }

                LayoutCategoryData category = getCategoryData(categoryId);
                if (category == null && categoryId.isEmpty()) {
                    category = retrieveBackupCategory();
                }

                if (category != null && !category.getLayouters().contains(algorithm)) {
                    category.getLayouters().add(algorithm);
                }
            }

            // Apply layout option dependencies registered with this registry (contrary to the code above, this code
            // requires that layout options we depend on have already been registered)
            for (Triple dep : optionDependencies) {
                LayoutOptionData source = getOptionData(dep.firstId);
                LayoutOptionData target = getOptionData(dep.secondId);
                if (source != null && target != null) {
                    source.getDependencies().add(Pair.of(target, dep.value));
                }
            }
            optionDependencies.clear();

            // Apply support information for supported layout options (contrary to the code above the code above, this
            // code requires that layout options we want to support have already been registered)
            for (Triple sup : optionSupport) {
                LayoutAlgorithmData algorithm = getAlgorithmData(sup.firstId);
                LayoutOptionData option = getOptionData(sup.secondId);
                if (algorithm != null && option != null) {
                    algorithm.addKnownOption(option, sup.value);
                }
            }
            optionSupport.clear();
        }

        /**
         * Returns the "Other" category. If there is none yet, creates one.
         *
         * @return the "Other" category.
         */
        private LayoutCategoryData retrieveBackupCategory() {
            LayoutCategoryData otherCategory = layoutCategoryMap.get("");
            if (otherCategory == null) {
                otherCategory = new LayoutCategoryData.Builder().id("").name("Other").create();
                layoutCategoryMap.put("", otherCategory);
            }

            return otherCategory;
        }

    }

}
