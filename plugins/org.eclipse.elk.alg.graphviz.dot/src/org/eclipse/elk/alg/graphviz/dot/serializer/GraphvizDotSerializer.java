/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.dot.serializer;

import java.io.IOException;
import java.util.List;

import org.eclipse.emf.ecore.EObject;
import org.eclipse.xtext.formatting2.FormatterPreferenceKeys;
import org.eclipse.xtext.formatting2.FormatterRequest;
import org.eclipse.xtext.formatting2.IFormatter2;
import org.eclipse.xtext.formatting2.regionaccess.ITextRegionAccess;
import org.eclipse.xtext.formatting2.regionaccess.ITextReplacement;
import org.eclipse.xtext.preferences.MapBasedPreferenceValues;
import org.eclipse.xtext.resource.SaveOptions;
import org.eclipse.xtext.serializer.impl.Serializer;

import com.google.common.collect.Maps;
import com.google.inject.Inject;
import com.google.inject.Provider;

/**
 * Specialization of the Xtext serializer that configures auto-linewrap after 255 characters.
 */
public class GraphvizDotSerializer extends Serializer {
    
    private static final int MAX_LINE_LENGTH = 255;
    
    @Inject
    private Provider<IFormatter2> formatter2Provider;

    @Inject
    private Provider<FormatterRequest> formatterRequestProvider;
    
    @Override
    protected void serialize(EObject obj, Appendable appendable, SaveOptions options) throws IOException {
        ITextRegionAccess regionAccess = serializeToRegions(obj);
        FormatterRequest request = formatterRequestProvider.get();
        MapBasedPreferenceValues preferences = new MapBasedPreferenceValues(Maps.<String, String> newLinkedHashMap());
        preferences.put(FormatterPreferenceKeys.maxLineWidth, MAX_LINE_LENGTH);
        request.setPreferences(preferences);
        request.setFormatUndefinedHiddenRegionsOnly(!options.isFormatting());
        request.setTextRegionAccess(regionAccess);
        IFormatter2 formatter2 = formatter2Provider.get();
        List<ITextReplacement> replacements = formatter2.format(request);
        regionAccess.getRewriter().renderToAppendable(replacements, appendable);
    }

}
