/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.grandom.validation

import org.eclipse.elk.core.debug.grandom.gRandom.Configuration
import org.eclipse.elk.core.debug.grandom.gRandom.Form
import org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage
import org.eclipse.xtext.validation.Check

/**
 * This class contains custom validation rules. 
 *
 * See https://www.eclipse.org/Xtext/documentation/303_runtime_concepts.html#validation
 */
class GRandomValidator extends AbstractGRandomValidator {
	
	static val ERRORS = newArrayList(
            Pair.of([Configuration c | c.form != Form.TREES && c.MD], 
                Pair.of("maxDegree only defined on trees.", GRandomPackage.Literals.CONFIGURATION__MD)),
            Pair.of([Configuration c | c.form != Form.TREES && c.MW], 
                Pair.of("maxWidth only defined on trees.", GRandomPackage.Literals.CONFIGURATION__MW)), 
            Pair.of([Configuration c | c.form != Form.BIPARTITE && c.PF],
                Pair.of("Partition Fraction only defined on bipartite graphs.", 
                    GRandomPackage.Literals.CONFIGURATION__PF
                ))
	);
	
	@Check
    def optionsRestrictions(Configuration conf){
        for (p : GRandomValidator.ERRORS) {
            if (p.key.apply(conf)) {
               error(p.value.key, p.value.value);
           }
        }
    }
	
}
