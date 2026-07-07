/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.ui;

import java.util.HashMap;
import java.util.Map;

import org.eclipse.elk.core.data.ILayoutMetaData;
import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutCategoryData;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.jface.viewers.ITreeContentProvider;
import org.eclipse.jface.viewers.Viewer;

/**
 * Content provider for displaying layout algorithms sorted by layout type.
 *
 * @author msp
 */
public class AlgorithmContentProvider implements ITreeContentProvider {

    /** the layout services used for this provider. */
    private LayoutMetaDataService layoutDataService;
    /** the filter map that stores visibility information. */
    private final Map<Object, Boolean> filterMap = new HashMap<Object, Boolean>();
    /** the current filter value. */
    private String filterValue;
    /** the current best filter match. */
    private ILayoutMetaData bestFilterMatch;
    
    @Override
    public void inputChanged(final Viewer viewer, final Object oldInput, final Object newInput) {
    }

    @Override
    public Object[] getElements(final Object inputElement) {
        if (inputElement instanceof LayoutMetaDataService) {
            layoutDataService = (LayoutMetaDataService) inputElement;
        }
        return layoutDataService.getCategoryData().toArray();
    }

    @Override
    public Object[] getChildren(final Object parentElement) {
        if (parentElement instanceof LayoutCategoryData) {
            LayoutCategoryData typeData = (LayoutCategoryData) parentElement;
            return typeData.getLayouters().toArray();
        }
        return null;
    }

    @Override
    public Object getParent(final Object element) {
        if (element instanceof LayoutAlgorithmData) {
            LayoutAlgorithmData layouterData = (LayoutAlgorithmData) element;
            return layoutDataService.getCategoryData(layouterData.getCategoryId());
        }
        return null;
    }

    @Override
    public boolean hasChildren(final Object element) {
        if (element instanceof LayoutCategoryData) {
            return ((LayoutCategoryData) element).getLayouters().size() > 0;
        }
        return false;
    }

    @Override
    public void dispose() {
    }
    
    /**
     * Update the filter value for this provider.
     * 
     * @param filter the new filter value
     */
    public void updateFilter(final String filter) {
        this.filterValue = filter;
        if (filterValue != null) {
            filterValue = filterValue.toLowerCase();
        }
        filterMap.clear();
        bestFilterMatch = null;
    }
    
    /**
     * Apply the current filter to the given element.
     * 
     * @param element an element from the content
     * @return true if the filter admits the element
     */
    public boolean applyFilter(final Object element) {
        Boolean result = filterMap.get(element);
        if (result == null) {
            if (filterValue != null && filterValue.length() > 0) {
                if (element instanceof LayoutCategoryData) {
                    LayoutCategoryData typeData = (LayoutCategoryData) element;
                    result = !typeData.getLayouters().isEmpty();
                    if (result) {
                        String typeName = typeData.getName().toLowerCase();
                        result = typeName.contains(filterValue);
                        if (result) {
                            for (LayoutAlgorithmData layouterData : typeData.getLayouters()) {
                                filterMap.put(layouterData, Boolean.TRUE);
                            }
                            if (bestFilterMatch == null
                                    || typeName.startsWith(filterValue)
                                    && !bestFilterMatch.getName().toLowerCase()
                                        .startsWith(filterValue)) {
                                bestFilterMatch = typeData;
                            }
                        } else {
                            boolean hasFilteredChild = false;
                            for (LayoutAlgorithmData layouterData : typeData.getLayouters()) {
                                hasFilteredChild |= applyFilter(layouterData);
                            }
                            result = hasFilteredChild;
                        }
                    }
                } else if (element instanceof LayoutAlgorithmData) {
                    LayoutAlgorithmData layouterData = (LayoutAlgorithmData) element;
                    String layouterName = layouterData.getName().toLowerCase();
                    if (layouterName.contains(filterValue)) {
                        result = Boolean.TRUE;
                        if (bestFilterMatch == null
                                || layouterName.startsWith(filterValue)
                                && !bestFilterMatch.getName().toLowerCase()
                                    .startsWith(filterValue)) {
                            bestFilterMatch = layouterData;                            
                        }
                    } else {
                        String bundle = layouterData.getBundleName();
                        result = bundle != null && bundle.toLowerCase().contains(filterValue);
                    }
                }
            } else if (element instanceof LayoutCategoryData) {
                result = !((LayoutCategoryData) element).getLayouters().isEmpty();
            } else {
                result = Boolean.TRUE;
            }
            filterMap.put(element, result);
        }
        return result;
    }
    
    /**
     * Returns the best match of the currently active filter.
     * 
     * @return the best filter match
     */
    public ILayoutMetaData getBestFilterMatch() {
        return bestFilterMatch;
    }

}
